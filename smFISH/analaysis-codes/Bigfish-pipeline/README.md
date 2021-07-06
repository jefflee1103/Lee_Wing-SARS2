# smFISH_Pypeline
Processing and analysis of smFISH data using Python repos (cellpose and bigfish)

## Requirements
Several python packages are required and need to be installed before running the script:
```python -m pip install yaml tifffile big-fish cellpose```
Do note that Cellpose will download extra files when it's ran for the first time. To avoid this, open a python interpreter shell and run:
```from cellpose import models```

## Running the script
The script is controlled with the `smFISH_analysis_config` YAML file. The intention is that the data is stored in an input directory specified by the `input_pattern` configuration parameter. All images recognised by the pattern will be processed. The output data will be stored in a new directory, located at the `output_dir` configuration parameter.