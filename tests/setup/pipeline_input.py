import unittest
from pathlib import Path
from src.setup.pipeline_input import FileFinder
from src.setup.config_loader import Config, StatisticalSettings, ReferencePaths, PipelineSettings, Paths

def create_base_config():
    return Config(
        version=1.0,
        run_name="test_run",
        paths=Paths(
            input_dir=Path("/fake/directory"),
            run_dir=Path("/fake/directory/run"),
            hg19=ReferencePaths(),
            hg38=ReferencePaths(),
        ),
        pipeline_settings=PipelineSettings(
            reference_genome="hg19",
            hicpro_raw_dirname="raw",
            hicpro_norm_dirname="iced",
            interaction_type="intra",
            iced_data=False,
            round_iced_matrices=False,
            intra_resolutions=[40000, 1000000, 20000],
            inter_resolutions=[],
            filter_blacklist=False,
            filter_cytobands=False,
            filter_self_interactions=False,
            remove_chromosomes=[],
            select_chromosomes=[],
            make_plots=False,
            select_specific_regions=False,
            omit_specific_regions=False,
            use_interaction_distance_filters=False,
            interaction_distance_filters={},
            output_format="csv",
            output_type="intra",
        ),
        statistical_settings=StatisticalSettings(
            spline_passes=1,
            fdr_threshold=0.01,
            metabin_occupancy=3,
            use_hicpro_bias=False,
            bias_lower_bound=0.1,
            bias_upper_bound=1.0,
            use_filtered_data_for_average_contact_probability=False,
            use_sequential_fdr=False,
        ),
    )

def create_iced_config():
    config = create_base_config()
    config.pipeline_settings.iced_data = True
    return config

def create_bias_config():
    config = create_base_config()
    config.statistical_settings.use_hicpro_bias = True
    return config

def create_iced_bias_config():
    config = create_iced_config()
    config.statistical_settings.use_hicpro_bias = True
    return config

class TestHicProInputFilePreparer(unittest.TestCase):

    def test_filter_files_on_resolution_raw(self):
        config = create_base_config()
        self.run_filter_test(config)

    def test_filter_files_on_resolution_iced(self):
        config = create_iced_config()
        self.run_filter_test(config)

    def test_filter_files_on_resolution_bias(self):
        config = create_bias_config()
        self.run_filter_test(config)

    def test_filter_files_on_resolution_iced_bias(self):
        config = create_iced_bias_config()
        self.run_filter_test(config)

    def run_filter_test(self, config):
        file_finder = FileFinder(config)

        test_cases = [
            ("MCF7_rep1_GSM4097072_40000.matrix", 40000),
            ("SampleX_1000000.matrix", 1000000),
            ("Sample_underscore_in_name_40000.matrix", 40000),
            ("Sample_with_100000_in_name_40000.matrix", 40000),
            ("Sample_100_with_12_lots_of_numbers_10_20000.matrix", 20000),
            ("hic_data/mcf10/matrix/mcf10/raw/1000000/mcf10_1000000.matrix", 1000000),
            ("hic_data/mcf10/matrix/mcf10/raw/1000000/mcf10_1000000_abs.bed", 1000000),  
            ("hic_DATA/mcF10/MatRiX/mcf10/raW/1000000/mCf10_1000000_abs.bed", 1000000),  
            ("hic_data/mcf10/matrix/mcf10/iced/1000000/mcf10_1000000_iced.matrix", 1000000),
            ("hic_data/mcf10/matrix/mcf10/iced/1000000/mcf10_1000000_iced.matrix.biases", 1000000),
        ]

        input_files = [Path(filename) for filename, _ in test_cases]
        
        filtered_files, found_resolutions = file_finder._filter_files_on_resolution(input_files)

        expected_resolutions_map = {Path(filename).name: resolution for filename, resolution in test_cases}

        extracted_resolutions = []
        for file in filtered_files:
            resolution = expected_resolutions_map[file.name]
            extracted_resolutions.append(resolution)
            print(f"File: {file.name}, Matched Resolution: {resolution}, Expected: {expected_resolutions_map[file.name]}")

        print("Extracted Resolutions:", extracted_resolutions)
        print("Expected Resolutions:", [resolution for _, resolution in test_cases])

        expected_resolutions_set = set(expected_resolutions_map.values())

        self.assertEqual(found_resolutions, expected_resolutions_set)


if __name__ == "__main__":
    unittest.main()