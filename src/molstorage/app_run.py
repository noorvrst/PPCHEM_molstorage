# src/molstorage/app.py

import os

class App:
    """A class to run the molstorage Streamlit interface."""

    def __init__(self) -> None:
        pass

    @staticmethod
    def run():
        """Starts the Streamlit interface."""
        dir_path = os.path.dirname(os.path.realpath(__file__))
        os.system(f"streamlit run {os.path.join(dir_path, 'app.py')}")