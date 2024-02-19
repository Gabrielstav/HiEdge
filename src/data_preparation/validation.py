# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import pandera as pa
from pandera.typing.dask import Series
from src.setup.config_loader import Config

class BedSchema(pa.DataFrameModel):
    chr: Series[str] = pa.Field(nullable=False)
    start: Series[int] = pa.Field(nullable=False)
    end: Series[int] = pa.Field(nullable=False)
    idx: Series[int] = pa.Field(nullable=False)

def get_matrix_schema(config):
    """
    Dynamically generates a schema for the matrix dataframe based on the config settings
    (not possible to dynamically set field types in based on config settings in pandera).
    """
    interaction_count_dtype = float if config.pipeline_settings.iced_data else int

    # Define the schema dynamically
    schema = pa.DataFrameSchema({
        "id_1": pa.Column(int),
        "id_2": pa.Column(int),
        "interaction_count": pa.Column(interaction_count_dtype),
    })

    return schema

def get_interaction_schema(config):
    """
    Also dynamically generate schema based on config for the interaction dataframe.
    :param config: Config instance
    :return: Interaction dataframe schema
    """

    interaction_count_dtype = float if config.pipeline_settings.iced_data else int

    # Define the schema dynamically
    schema = pa.DataFrameSchema({
        "chr_1": pa.Column(str),
        "start_1": pa.Column(int),
        "end_1": pa.Column(int),
        "chr_2": pa.Column(str),
        "start_2": pa.Column(int),
        "end_2": pa.Column(int),
        "interaction_count": pa.Column(interaction_count_dtype),
        "idx_1": pa.Column(int),
        "idx_2": pa.Column(int)
    })

    return schema


def validate_bed_data(df):
    validated_df = BedSchema.validate(df, lazy=True)
    return validated_df

def validate_matrix_data(df, config: Config):
    schema = get_matrix_schema(config)
    return schema.validate(df, lazy=True)

def validate_interaction_data(df, config: Config):
    schema = get_interaction_schema(config)
    return schema.validate(df, lazy=True)
