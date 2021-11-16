import os
from collections import namedtuple, OrderedDict


class Pipeline:
    """
    Stores information of samples, their references and targeted panels in structures suitable
    for SnakeMake expand() function
    """
    WGS_PANEL = 'wgs'
    
    def __init__(self):
        self.samples = []
        self.metadata = []
        self.references = []
        self.panels = []
        self.sample_references = []
        self.__reference_map = {}
    
    def add(self, samples, metadata, reference, panel, prebuilt_reference, would_be_downloaded):
        """
        Add samples intended for analysis.
        :param samples: list of sample names to be analysed
        :param metadata: list of sample descriptions
        :param reference: name of reference fasta file for samples
        :param panel: name of targeted panel for samples
        :param prebuilt_reference: if reference is only a parameter to a specific tool that utilize its own database,
                                   and so reference/{reference}/{reference}.fa could not exist yet
        :param would_be_downloaded: if reference would be created in the Snakelines execution,
                                   and so reference/{reference}/{reference}.fa could not exist yet
        """
        
        # Check, if reference file exists
        if reference and not (prebuilt_reference or would_be_downloaded):
            fasta = 'reference/{reference}/{reference}.fa'.format(reference=reference)
            assert os.path.exists(fasta), 'Reference fasta {} does not exist'.format(fasta)
        
        panel_str = self._panel_str(panel, reference)
        
        if panel_str not in self.panels:
            self.panels.append(panel_str)
        
        # Extend lists
        self.samples.extend(samples)
        self.metadata.extend(metadata)

        if reference:
            if reference not in self.references:
                self.references.append(reference)
            
            if reference in self.__reference_map:
                self.__reference_map[reference].extend(samples)
            else:
                self.__reference_map[reference] = samples
            
            for sample in samples:
                # TODO here should be check if tuple already stored
                SampleReference = namedtuple('SampleReference', 'sample reference panel')
                self.sample_references.append(SampleReference(sample, reference, panel_str))
    
    def _panel_str(self, panel: OrderedDict, reference: str) -> str:
        # select panel
        if panel and panel['name'] != self.WGS_PANEL:
            # defined panel requires a bed file
            assert reference, 'Panel cannot be defined without a reference'
            
            # define original bed
            panel_str = panel['name']
            bed = 'reference/{reference}/annotation/{panel}/regions.bed' \
                .format(reference=reference, panel=panel_str)
            
            if 'flank' in panel:
                # use flanked bed
                panel_str = '{panel}_flank{flank}'.format(panel=panel['name'], flank=panel['flank'])
                flanked_bed = 'reference/{reference}/annotation/{panel}/regions.bed' \
                    .format(reference=reference, panel=panel_str)
                
                if not os.path.exists(bed) and not os.path.exists(flanked_bed):
                    # neither flanked bed or original bed exist
                    raise FileNotFoundError('Panel bed file %s does not exist' % bed)
            
            elif not os.path.exists(bed):
                raise FileNotFoundError('Panel bed file %s does not exist' % bed)
        
        else:
            # undefined panel is whole genome "panel"
            panel_str = self.WGS_PANEL
            
        return panel_str
    
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
