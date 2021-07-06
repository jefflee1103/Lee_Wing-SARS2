import pathlib
import glob
import signal
import queue
import multiprocessing
import threading
import time
import getpass
import omero.gateway
import yaml
import scipy
from skimage.morphology import white_tophat, black_tophat, disk
import numpy as np
import tifffile
import bigfish.stack as stack
import bigfish.detection as detection
from cellpose import models


# calculate psf (thank you MK), with edit for consistent nomenclature
def calculate_psf(voxel_size_z, voxel_size_yx, Ex, Em, NA, RI, microscope):
    """
    Use the formula implemented in Matlab version (sigma_PSF_BoZhang_v1)
    to calculate the theoretical PSF size.
    """
    if microscope == "widefield":
        psf_yx = 0.225 * Em / NA
        psf_z = 0.78 * RI * Em / (NA ** 2)
    elif microscope in ("confocal", "nipkow"):
        psf_yx = 0.225 / NA * Ex * Em / np.sqrt(Ex ** 2 + Em ** 2)
        psf_z = 0.78 * RI / NA ** 2 * Ex * Em / np.sqrt(Ex ** 2 + Em ** 2)
    else:
        # Unrecognised microscope
        raise Exception(
            "Unrecognised microscope argument for function calculate_psf()."
        )
    return psf_z, psf_yx


# subtract background
def subtract_background(image, radius, light_bg=False):
    # you can also use 'ball' here to get a slightly smoother result at the
    # cost of increased computing time
    str_el = disk(radius)
    if light_bg:
        return black_tophat(image, str_el)
    else:
        return white_tophat(image, str_el)


def image_processing_function(image_loc, config):
    # Read the image into a numpy array of format ZCYX
    if isinstance(image_loc, str):
        image_name = pathlib.Path(image_loc).stem
        image = tifffile.imread(image_loc)
    else:
        # Establish connection with OMERO and actually connect
        conn = omero.gateway.BlitzGateway(
            host="omero1.bioch.ox.ac.uk",
            port=4064,
            # group=config["OMERO_group"],
            username=config["OMERO_user"],
            passwd=config["password"],
        )
        conn.connect()
        conn.SERVICE_OPTS.setOmeroGroup(-1)
        # Create a thread to keep the connection alive
        ka_thread = threading.Thread(
            target=keep_connection_alive, args=(conn,)
        )
        ka_thread.daemon = True
        ka_thread.start()
        # Derive image and its name
        image_name = pathlib.Path(image_loc[0]).stem
        remote_image = conn.getObject("Image", image_loc[1])
        image = np.array(
            list(
                remote_image.getPrimaryPixels().getPlanes(
                    [
                        (z, c, 0)
                        for z in range(0, remote_image.getSizeZ())
                        for c in range(0, remote_image.getSizeC())
                    ]
                )
            )
        )
        image = image.reshape(
            image.shape[0] // remote_image.getSizeC(),
            remote_image.getSizeC(),
            image.shape[1],
            image.shape[2],
        )

    # segment with cellpose
    seg_img = np.max(image[:, config["seg_ch"], :, :], 0)
    if config["cp_search_string"] in image_name:
        seg_img = np.clip(seg_img, 0, config["cp_clip"])
    seg_img = scipy.ndimage.median_filter(
        seg_img, size=config["median_filter"]
    )
    model = models.Cellpose(gpu=config["gpu"], model_type="cyto")
    channels = [0, 0]  # greyscale segmentation
    masks = model.eval(
        seg_img,
        channels=channels,
        diameter=config["diameter"],
        do_3D=config["do_3D"],
        flow_threshold=config["flow_threshold"],
        cellprob_threshold=config["cellprob_threshold"],
    )[0]

    # Calculate PSF
    psf_z, psf_yx = calculate_psf(
        config["voxel_size_z"],
        config["voxel_size_yx"],
        config["ex"],
        config["em"],
        config["NA"],
        config["RI"],
        config["microscope"],
    )
    sigma = detection.get_sigma(
        config["voxel_size_z"], config["voxel_size_yx"], psf_z, psf_yx
    )

    for image_channel in config["channels"]:
        # detect spots
        rna = image[:, image_channel, :, :]
        # subtract background
        rna_no_bg = []
        for z in rna:
            z_no_bg = subtract_background(z, config["bg_radius"])
            rna_no_bg.append(z_no_bg)
        rna = np.array(rna_no_bg)

        # LoG filter
        rna_log = stack.log_filter(rna, sigma)

        # local maximum detection
        mask = detection.local_maximum_detection(rna_log, min_distance=sigma)

        # tresholding
        if image_channel == config["smFISH_ch1"]:
            threshold = config["smFISH_ch1_thresh"]
        elif image_channel == config["smFISH_ch2"]:
            threshold = config["smFISH_ch2_thresh"]
        else:
            print("smFISH channel and threshold not correctly defined!")

        spots, _ = detection.spots_thresholding(rna_log, mask, threshold)

        # detect and decompose clusters
        spots_post_decomposition = detection.decompose_cluster(
            rna,
            spots,
            config["voxel_size_z"],
            config["voxel_size_yx"],
            psf_z,
            psf_yx,
            alpha=config["alpha"],  # impacts number of spots per cluster
            beta=config["beta"],  # impacts the number of detected clusters
        )[0]

        # separate spots from clusters
        spots_post_clustering, foci = detection.detect_foci(
            spots_post_decomposition,
            config["voxel_size_z"],
            config["voxel_size_yx"],
            config["bf_radius"],
            config["nb_min_spots"],
        )

        # extract cell level results
        image_contrasted = stack.rescale(rna, channel_to_stretch=0)
        image_contrasted = stack.maximum_projection(image_contrasted)
        rna_mip = stack.maximum_projection(rna)

        fov_results = stack.extract_cell(
            cell_label=masks.astype(np.int64),
            ndim=3,
            rna_coord=spots_post_clustering,
            others_coord={"foci": foci},
            image=image_contrasted,
            others_image={"smfish": rna_mip},
        )

        # save bigfish results
        for i, cell_results in enumerate(fov_results):
            output_path = pathlib.Path(config["output_dir"]).joinpath(
                f"{image_name}_ch{image_channel + 1}_results_cell_{i}.npz"
            )
            stack.save_cell_extracted(cell_results, str(output_path))

        # save reference spot for each image
        # (Using undenoised image! not from denoised!)
        reference_spot_undenoised = detection.build_reference_spot(
            rna,
            spots,
            config["voxel_size_z"],
            config["voxel_size_yx"],
            psf_z,
            psf_yx,
            alpha=config["alpha"],
        )
        spot_output_path = pathlib.Path(config["output_refspot_dir"]).joinpath(
            f"{image_name}_reference_spot_ch{image_channel + 1}"
        )
        stack.save_image(
            reference_spot_undenoised, str(spot_output_path), "tif"
        )

    # Close the OMERO connection
    conn.close()


def worker_function(jobs, results):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while not jobs.empty():
        try:
            job = jobs.get(block=False)
            results.put(image_processing_function(*job))
        except queue.Empty:
            pass


def keep_connection_alive(conn):
    while True:
        conn.keepAlive()
        time.sleep(60)


def main():
    jobs = multiprocessing.Queue()
    results = multiprocessing.Queue()

    # Load the config file
    with open("smFISH_analysis_config.yaml") as fi:
        config = yaml.load(fi, Loader=yaml.Loader)

    # Check if output directories exists; try to create them if they don't
    pathlib.Path(config["output_dir"]).mkdir(exist_ok=True)
    pathlib.Path(config["output_refspot_dir"]).mkdir(exist_ok=True)

    # Fill the job queue either with local or remote files
    if "OMERO_user" in config:
        # Ask for password
        password = getpass.getpass(
            f"Type password for user '{config['OMERO_user']}':"
        )
        config["password"] = password
        # Establish connection with OMERO and actually connect
        conn = omero.gateway.BlitzGateway(
            host="omero1.bioch.ox.ac.uk",
            port=4064,
            # group=config["OMERO_group"],
            username=config["OMERO_user"],
            passwd=password,
        )
        conn.connect()
        conn.SERVICE_OPTS.setOmeroGroup(-1)
        # Fetch the images and their original names
        for dataset_id in config["OMERO_datasets"]:
            for image in conn.getObject("dataset", dataset_id).listChildren():
                for orig_file in image.getImportedImageFiles():
                    jobs.put(((orig_file.getName(), image.getId()), config))
    else:
        # Get images using the input path pattern
        image_paths = glob.glob(config["input_pattern"])
        for image_path in image_paths:
            jobs.put((image_path, config))

    # Close OMERO connection
    conn.close()

    # Start workers
    workers = []
    for i in range(config["number_of_workers"]):
        p = multiprocessing.Process(
            target=worker_function, args=(jobs, results)
        )
        p.start()
        workers.append(p)

    # Wait for workers to complete
    try:
        for worker in workers:
            worker.join()
    except KeyboardInterrupt:
        for worker in workers:
            worker.terminate()
            worker.join()


if __name__ == "__main__":
    # Set the process start method; spawn is default on Windows
    multiprocessing.set_start_method("spawn")
    # Call the main function
    main()
