from typing import Union
import os

def get_value_from_dict(vocab: dict, key: str) -> Union[str, list, int]:
    try:
        value = vocab[key]
        if value is None:
            raise ValueError(f'Value for key {key} is None in row {vocab}')
        return value
    except KeyError:
        raise KeyError(f"Row {vocab} does not have a key {key} to get value")
    
def get_files_with_extension(directory, extension):
    files_with_extension = []
    for filename in os.listdir(directory):
        if filename.endswith(extension):
            files_with_extension.append(os.path.join(directory, filename))
    return files_with_extension