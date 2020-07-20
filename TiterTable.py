#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 00:09:51 2020

@author: avicenna
"""

import numpy as np
import pandas as pd
import AntigenicDatabase as ad
from tools import MAPTYPES, convert_titer

    
class TiterTable():
    
    '''
    Titer table object which can be initialized from results, serum, antigen 
    json files. Inside each titer table an AntigenDataset and SerumDataset
    of the antigens and sera of the titer table is initialized.
    
    During initialization convert_to_dict creates 4 dictionaries 
    numerical_titer_dict, titer_dict, lessthan_dict, greaterthan_dict
    that respectively map (antigen_id, serum_id) pairs to numerical titer values, 
    titer values as it appears in the results database and to bool values that 
    tell whether or  not the titer is thresholded. 
    
    At the cost of reduntant serum_id and antigen_id information stored in 
    keys, this grants precision and allows one to manipulate data securely
    without worrying about messing up order of sera or antigens. For instance
    if you want to create a table with antigens and sera ordered in a specific way, 
    all you have to do is get two arrays: antigen_ids, serum_ids and loop through
    these ids while retrieving the relevant titer from dictionaries and recording
    it to your table.
    
    '''
    
    def __init__(self, data, serum_json, antigen_json):
        
        assert isinstance(data, MAPTYPES), 'Provided data should be a map like type'

        self.data = data
        if 'file' in data:
            self.file_name = data['file']
        self.serum_ids = data['serum_ids']
        self.antigen_ids = data['antigen_ids']
        
        self.serum_dataset = ad.SerumDataset(serum_json, id_subset = self.serum_ids)
        self.antigen_dataset = ad.AntigenDataset(antigen_json, id_subset = self.antigen_ids)
        
        
        self.antigens = self.antigen_dataset.entries
        self.sera =self.serum_dataset.entries
        
        self.antigen_ids = [x.id for x in self.antigens]
        self.serum_ids = [x.id for x in self.sera]
        self.serum_longs = [x.long for x in self.sera]
        self.antigen_longs = [x.long for x in self.antigens]
        self.serum_strain_ids = [x.strain_id for x in self.sera]
        
        self.antigen_shorts = []
        self.sera_shorts = []
        
        self.titers = data['titers']

        self.convert_to_dicts()
        
        self.general_health_check()
        
        if len(self.serum_strain_ids) != len(set(self.serum_strain_ids)):
            print('There may be repeated measurements in this dataset (repeated serum strain ids).')
   
    def general_health_check(self):
    
        number_of_antigens, number_of_sera = np.array(self.titers).shape
        
        assert (number_of_antigens, number_of_sera) == (len(self.antigen_ids), len(self.serum_ids)), 'Table size {} should be consistent with number of antigens {} and sera ids {}'.format((number_of_antigens, number_of_sera), len(self.antigen_ids), len(self.serum_ids))                                        
        assert (number_of_antigens, number_of_sera) == (len(self.antigen_longs), len(self.serum_longs)), 'Table size should be consistent with number of antigens and sera longs'
        
        
        assert (number_of_antigens == len(self.titer_dict) 
                and len(self.titer_dict) == len(self.numerical_titer_dict)
                and len(self.titer_dict) == len(self.greaterthan_dict)
                and len(self.titer_dict) == len(self.lessthan_dict)), 'All dictionary elements should have the same length which is equal to number of antigens'
        
        td_length = list(set([len(self.titer_dict[x]) for x in self.titer_dict]))
        ntd_length = list(set([len(self.numerical_titer_dict[x]) for x in self.numerical_titer_dict]))
        gtd_length = list(set([len(self.greaterthan_dict[x]) for x in self.numerical_titer_dict]))
        ltd_length = list(set([len(self.lessthan_dict[x]) for x in self.numerical_titer_dict]))
        
        assert len(td_length) == 1, 'Each element of titer dict should have the same length'
        assert len(ntd_length) == 1, 'Each element of numerical titer dict should have the same length'
        assert len(gtd_length) == 1, 'Each element of greaterthans dict should have the same length'
        assert len(ltd_length) == 1, 'Each element of lessthans titer dict should have the same length'
        
        assert (number_of_sera == td_length[0] 
                and td_length[0] == ntd_length[0] 
                and ntd_length[0] == gtd_length[0] 
                and gtd_length[0] == ltd_length[0]), 'All dictionaries elements should contain same number of elements which should equal the number of sera'
                        

        
        
    def convert_to_dicts(self):
        
        self.numerical_titer_dict = {}
        self.lessthan_dict = {}
        self.greaterthan_dict = {}
        self.titer_dict = {}
        
        for antigen_id, titer_entry in zip(self.antigen_ids, self.titers):
            
            self.lessthan_dict[antigen_id] = { sid:(True if isinstance(x,str) and '<' in x else False) for sid,x in zip(self.serum_ids, titer_entry) }
            self.greaterthan_dict[antigen_id] = { sid:(True if isinstance(x,str) and '>' in x else False) for sid,x in zip(self.serum_ids, titer_entry) }
            
            self.numerical_titer_dict[antigen_id] = {sid:convert_titer(x) for sid,x in zip(self.serum_ids, titer_entry) }
            self.titer_dict[antigen_id] = {sid:x for sid,x in zip(self.serum_ids, titer_entry)}
    
    def homologous_sera_order(self):
            
        homologous_order = []
        
        for antigen_id in self.antigen_ids:
            
            homologous_order.append(self.serum_strain_ids.index(antigen_id))
    
    def to_df(self, as_is=False, do_rounding=False, thresholded=False,
              serum_order_ids=None, antigen_order_ids=None, extra_rows = None,
              extra_columns = None, add_ids = False, add_serum_strain_ids = False,
              antigen_names = None, serum_names = None):
        
        ''' 
        Produces a dataframe from titer dictionaries
        
        Parameters
        -----------      
        as_is: If True, the values are written including /, < and >
        This overrides all the parameters below.
        
        do_rounding: If True, non integer values are rounded to integers
        
        thresholded: If True, thresholded values are recorded as <,> if not
        < is recorded as 1/2 times the value and > is recorded as 2 times the value
        
        serum_order_ids and antigen_order_ids should be list of antigen and sera ids from 
        the TiterTable. Then the dataframe will be ordered according to these.
        
        extra_rows and extra_columns are dictionaries which are added to the table
        at the beginning and should have the format {name1:data1, ...}
        
        If add_ids and add_serum_strain_ids are True, then these are added as extra rows
        and columns (serum_strain_ids only as rows )
        
        If antigen_names and serum_names are provided these are used as index
        and column names for the dataframe
        
        '''      
        
        if antigen_order_ids is None:
            antigen_ids = self.antigen_ids
        else:
            assert set(antigen_order_ids).issubset(set(self.antigen_ids)),'Antigen order ids should be a subset of table antigen ids'
            antigen_ids = antigen_order_ids
        
        if serum_order_ids is None:
            serum_ids = self.serum_ids
        else:
            assert set(serum_order_ids).issubset(set(self.serum_ids)),'Serum order ids should be a subset of table serum ids'
            serum_ids = serum_order_ids
        
        
        if extra_rows is None:
            extra_rows = {}
        else:
            assert isinstance(extra_rows, MAPTYPES)
            assert all( len(extra_rows[x])==len(serum_ids) for x in extra_rows), 'extra rows should have the same length as number of sera'
                  
            
        if extra_columns is None:
            extra_columns = {}
        else:
            assert isinstance(extra_columns, MAPTYPES)
            assert all( len(extra_columns[x])==len(antigen_ids) for x in extra_columns), 'extra columnss should have the same length as number of antigens'
            
            
        if add_ids:
            extra_columns['id'] = antigen_ids
            extra_rows['id'] = serum_ids
            
        if add_serum_strain_ids:
            extra_rows['serum strain id'] = [self.serum_dataset.get_entry(x)[0].strain_id for x in serum_ids]
            
            
        if antigen_names is None:
            antigen_names =  [self.antigen_dataset.get_entry(x)[0].long for x in antigen_ids ] 
        if serum_names is None:
            serum_names = [self.serum_dataset.get_entry(x)[0].long for x in serum_ids ] 
            
        columns = list(extra_columns.keys()) + serum_names
        indices = list(extra_rows.keys()) + antigen_names
        df = pd.DataFrame(columns=columns, index=indices)
        raw_data = pd.DataFrame(columns=serum_names, index=antigen_names)
        
        
        
        c_start = len(extra_columns)
        r_start = len(extra_rows)
        
        
        
        for ind,key in enumerate(extra_columns):
            df.iloc[r_start:,ind] = extra_columns[key]
        
        for ind,key in enumerate(extra_rows):
            df.iloc[ind, c_start:] = extra_rows[key]
        
        for col_i,serum_id in enumerate(serum_ids):
            for row_i,antigen_id in enumerate(antigen_ids):
                
                if not as_is:
                    val = self.numerical_titer_dict[antigen_id][serum_id]
                    
                    if do_rounding:
                        val = int(np.round(val))
                          
                    if thresholded:
                        if self.lessthan_dict[antigen_id][serum_id]:
                            val = '<' + str(2 * val)
                            
                        elif self.greaterthan_dict[antigen_id][serum_id]:
                            val = '>' + str(val / 2)
                            
                    else:
                        
                        val = str(val)
                else:
                    val = self.titer_dict[antigen_id][serum_id]
        
                df.iloc[row_i + r_start, col_i + c_start] = val
                raw_data.iloc[row_i, col_i] = val
        
        
    
        return df, raw_data
        
    
    def __str__(self):
        
        return 'Titer data extracted from {}: '.format(self.file_name) 

