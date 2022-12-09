from abc import ABC, abstractmethod
from typing import List

from burden.burden_ingester import BurdenAssociationPack


class ToolRunner(ABC):

    def __init__(self, association_pack: BurdenAssociationPack, output_prefix: str):
        self._association_pack = association_pack
        self._output_prefix = output_prefix
        self._outputs = []

    def get_outputs(self) -> List[str]:
        return self._outputs

    @abstractmethod
    def run_tool(self) -> None:
        pass
