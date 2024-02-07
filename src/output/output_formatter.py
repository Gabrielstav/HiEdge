# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from pathlib import Path


class OutputConfigurator:
    def __init__(self, config):
        self.config = config

    def run(self, input_data):
        output_type = self.config.output_type
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
        columns = ["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "q_value", "confidence_estimate", "interaction_count"]
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
        file_name = f"{input_data.metadata.experiment}_{input_data.metadata.resolution}_{suffix}.csv"
        output_dir = Path(self.config.output_dir) / interaction_type
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir / file_name

    def _write_to_file(self, df, file_path):
        output_format = self.config.pipeline_settings.output_format
        file_name = file_path.name

        if output_format == "csv":
            df.to_csv(file_path, single_file=True)
        elif output_format == "parquet":
            df.to_parquet(file_path)
        elif output_format == "hdf5":
            # HDF5 writing needs special handling for Dask
            df.compute().to_hdf(file_path, key='data', mode='w')
        elif output_format == "txt":
            df.to_csv(file_path, single_file=True, sep='\t', header=False, index=False)
        else:
            raise ValueError(f"Unsupported output format: {output_format}")

        print(f"Output file {file_name} written to {file_path}")



