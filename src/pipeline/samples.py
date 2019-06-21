import os
from collections import namedtuple


class Pipeline:
    """
    Stores information of samples, their references and targeted panels in structures suitable
    for SnakeMake expand() function
    """

    def __init__(self):
        self.samples = []
        self.references = []
        self.panels = []
        self.sample_references = []
        self.sample_defs = []
        self.__reference_map = {}

    def add(self, samples, reference, panel, prebuilt_reference):
        """
        Add samples intended for analysis.
        :param samples: list of sample names to be analysed
        :param reference: name of reference fasta file for samples
        :param panel: name of targeted panel for samples
        :param prebuilt_reference: if reference is only a parameter to a specific tool that utilize its own database,
                                   and so reference/{reference}/{reference}.fa does not exist
        """

        # Check, if reference file exists
        if reference and not prebuilt_reference:
            fasta = 'reference/{reference}/{reference}.fa'.format(reference=reference)
            assert os.path.exists(fasta), 'Reference fasta {} does not exist'.format(fasta)

        # Check, if panel bed file exists
        if panel:
            assert reference, 'Panel cannot be defined without reference'
            bed = 'references/{reference}/{reference}/annotation/{panel}/regions.bed'.format(reference=reference, panel=panel)
            assert os.path.exists(bed), 'Panel bed file {} does not exist'.format(bed)

        # Extend lists
        self.samples.extend(samples)

        if reference:

            if reference not in self.references:
                self.references.append(reference)

            if reference in self.__reference_map:
                self.__reference_map[reference].extend(samples)
            else:
                self.__reference_map[reference] = samples

            for sample in samples:
                # TODO here should be check if tuple already stored
                self.sample_references.append(namedtuple('SampleReference', 'sample reference')(sample, reference))

        if panel:
            if panel not in self.panels:
                self.panels.append(panel)
            for sample in samples:
                # TODO here should be check if tuple already stored
                self.sample_defs.append(namedtuple('SampleDef', 'sample reference panel')(sample, reference, panel))

    def is_empty(self):
        """
        Has at least one sample to be analysed?
        """
        return len(self.samples) == 0

    def samples_for(self, reference):
        """
        Samples that should be analysed with defined reference
        :param reference: name of the reference
        """
        return self.__reference_map[reference]

    def __str__(self):
        """
        Returns the string representation of this class
        :return: str - string representation
        """
        res = "Samples: " + ', '.join(self.samples) + '\n'
        res += "References: " + ', '.join(self.references) + '\n'
        res += "Panels: " + ', '.join(self.panels) + '\n'
        return res