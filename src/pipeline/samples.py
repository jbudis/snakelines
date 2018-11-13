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

    def add(self, samples, reference, panel):
        """
        Add samples intended for analysis.
        :param samples: list of sample names to be analysed
        :param reference: name of reference fasta file for samples
        :param panel: name of targeted panel for samples
        """

        # Check if both R1 and R2 files for all samples exist
        for sample in samples:
            r1, r2 = 'reads/original/{}_R1.fastq.gz'.format(sample), 'reads/original/{}_R2.fastq.gz'.format(sample)
            assert os.path.exists(r1), 'Read file {} does not exists'.format(r1)
            assert os.path.exists(r2), 'Read file {} does not exists'.format(r2)

        # Check, if reference file exists
        if reference:
            fasta = 'reference/{reference}/{reference}.fa'.format(reference=reference)
            assert os.path.exists(fasta), 'Reference fasta {} does not exists'.format(fasta)

        # Check, if panel bed file exists
        if panel:
            assert reference, 'Panel cannot be defined without reference'
            bed = 'references/{reference}/{reference}/annotation/{panel}/regions.bed'.format(reference=reference, panel=panel)
            assert os.path.exists(bed), 'Panel bed file {} does not exists'.format(bed)

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
