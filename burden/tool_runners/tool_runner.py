from abc import ABC, abstractmethod
from pathlib import Path
from typing import List

from burden.burden_ingester import BurdenAssociationPack
from general_utilities.mrc_logger import MRCLogger


class ToolRunner(ABC):

    def __init__(self, association_pack: BurdenAssociationPack, output_prefix: str):

        self._logger = MRCLogger(__name__).get_logger()
        self._association_pack = association_pack
        self._output_prefix = output_prefix
        self._outputs = []

    def get_outputs(self) -> List[Path]:
        return self._outputs

    @abstractmethod
    def run_tool(self) -> None:
        pass
