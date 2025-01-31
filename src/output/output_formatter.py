# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from pathlib import Path
from dask import dataframe as dd

class OutputConfiguratorRunner:
    def __init__(self, config, input_data):
        self.config = config
        self.input_data = input_data

    def run(self):
        output_configurator = OutputConfigurator(config=self.config)
        print(f"Running output configurator with input data: {self.input_data}")
        output_configurator.run(self.input_data)


class OutputConfigurator:
    def __init__(self, config):
        self.config = config

    def run(self, input_data):
        print(f"Head of input data in output configuration: {input_data.data.head(5)}")
        output_type = self.config.pipeline_settings.output_type
        if output_type == "default":
            self._write_default_output(input_data)
        elif output_type == "verbose":
            self._write_verbose_output(input_data)
        elif output_type == "qval":
            self._write_qval_output(input_data)
        else:
            raise ValueError(f"Unknown output type: {output_type}")

    def _write_default_output(self, input_data):
        file_path = self._construct_file_path(input_data, suffix="default")
        columns = ["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count", "p_value", "q_value"]
        self._write_to_file(input_data.data[columns], file_path)

    def _write_verbose_output(self, input_data):
        file_path = self._construct_file_path(input_data, suffix="verbose")
        self._write_to_file(input_data.data, file_path)

    def _write_qval_output(self, input_data):
        file_path = self._construct_file_path(input_data, suffix="qval")
        columns = ["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "q_value"]
        self._write_to_file(input_data.data[columns], file_path)

    def _construct_file_path(self, input_data, suffix):
        interaction_type = input_data.metadata.interaction_type

        if self.config.pipeline_settings.output_format == "csv":
            file_name = f"{input_data.metadata.experiment}_{input_data.metadata.resolution}_{suffix}.csv"
        elif self.config.pipeline_settings.output_format == "parquet":
            file_name = f"{input_data.metadata.experiment}_{input_data.metadata.resolution}_{suffix}.parquet"
        elif self.config.pipeline_settings.output_format == "hdf5":
            file_name = f"{input_data.metadata.experiment}_{input_data.metadata.resolution}_{suffix}.h5"
        elif self.config.pipeline_settings.output_format == "txt":
            file_name = f"{input_data.metadata.experiment}_{input_data.metadata.resolution}_{suffix}.txt"
        else:
            raise ValueError(f"Unsupported output format: {self.config.pipeline_settings.output_format}. Supported formats are csv, parquet, hdf5, txt.")

        output_dir = Path(self.config.paths.output_dir) / interaction_type
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir / file_name

    def _write_to_file(self, df, file_path):
        output_format = self.config.pipeline_settings.output_format
        file_name = file_path.name

        if not isinstance(df, dd.DataFrame):
            print(f"Output dataframe is not a Dask dataframe: {type(df)}")

        if output_format == "csv":
            df.to_csv(file_path, index=False)  # Added index=False here
        elif output_format == "parquet":
            df.to_parquet(file_path)
        elif output_format == "hdf5":
            # HDF5 writing needs special handling for Dask
            df.compute().to_hdf(file_path, key="data", mode="w")
        else:  # txt format
            df.to_csv(file_path, sep="\t", header=True, index=False)

        print(f"Output file {file_name} written to {file_path}")



