#!/usr/bin/env python

""" MultiQC module to parse output from MetaData """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import json
import os
import sys
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Metadata', anchor='metadata',
        info="metadata from seQc repos.")

        # Parse logs
        self.metadata = dict()
        for f in self.find_log_files('metadata', filehandles=True):
            self.parse_metadata(f)

        # Filter to strip out ignored sample names
        self.metadata = self.ignore_samples(self.metadata)

        if len(self.metadata) == 0:
            raise UserWarning

        log.info("Found {} logs".format(len(self.metadata)))
        self.write_data_file(self.metadata, 'multiqc_metadata')

        # Add drop rate to the general stats table
        headers = OrderedDict()
        headers['Project_ID'] = {
            'title': 'Project_ID',
            'description': 'Project_ID'
        }
        """
	headers['volume'] = {
            'title': 'Sample Volume',
            'description': 'Sample Volume',
            'max': 175,
            'min': 0
        }
	headers['conc'] = {
	    'title': 'Sample Concentration',
	    'description': 'Sample Concentration',
	    'max': 188.7,
	    'min': 0
	}
	headers['pool_id'] = {
	    'title': 'Pool ID',
	    'description': 'Mock seqencing pool'
	}
	headers['RQS'] = {
	    'title': 'RQS',
	    'description': 'RNA-Quality Score',
	    'max': 7,
	    'min': 0
	}
        """
        self.general_stats_addcols(self.metadata, headers)

        # Make barplot
        # self.metadata_boxplot()

    def parse_metadata(self, f):
         for line in f['f']:
             data = json.loads(line)
             s_name = os.path.abspath(f['f'].name).split('/')[-2]
             suffices = ['_R1', '_R2']
             for s in suffices:
                  self.metadata[s_name + s] = dict()
                  print(data)
                  self.metadata[s_name + s]['Project_ID'] = data["Project_ID"]
                  """
                  self.metadata[s_name + s]['volume'] = data["volume"]
                  self.metadata[s_name + s]['conc'] = data["conc"]
                  self.metadata[s_name + s]['pool_id'] = data['pool_id']
                  self.metadata[s_name + s]['RQS'] = data['RQS']
                  """
