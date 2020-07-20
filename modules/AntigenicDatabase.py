#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:33:08 2020

@author: Sina

A module for parsing antigen and serum datasets from acorg. 


"""
import numpy as np
import re
import datetime
import collections
import tools as t
import editdistance as ed
from itertools import compress
from CityList import CityList

city_list = CityList()
city_names = [x.name for x in city_list.cities]
now = datetime.datetime.now()
max_year=now.year
min_year=1800

MAPTYPES = (dict, collections.abc.Mapping)


class Entry:

    '''parent class from which SerumEntry and AntigenEntry classes 
    are constructed. This is either an entry about a serum or an antigen, 
    initialized from data which is of the form dictionary and should atleast
    containts id and long as keys. Supports deep dictionary search and aliased
    searched for matching queries with regular expressions. It also finds a 
    short name for the entry by abbreviating the city name.'''
    
    def __init__(self, data):
        
        assert isinstance(data, MAPTYPES), 'Provided data should be a map like type'

        self._data = data  # keep a copy of the data inside the class for deep searches, tests etc

        # read id, long name which are required to exist
        assert 'id' in self._data, 'Data should contain the key: id'
        self.id = data['id']
        assert 'long' in self._data, 'Data should contain the key: long'
        self.long = data['long']
        
        # store other mappings in the variable properties for easy access
        self.properties ={key:self._data[key] for key in self._data if key not in ['id', 'long']}
        self.short = self._find_short_name()



    def _find_short_name(self):
        
        short_name = self.long
        cities = list(set(city_list.search(short_name, split=True)))
        
        
        for city in cities:
            if city.name in short_name:
                short_name = short_name.replace(city.name, city.abb)
            else:
                split_name = short_name.split('/')
                edists = []
                for name in split_name:
                    edists.append(ed.eval(name, city.name))
                    
                split_name[np.argmin(np.array(edists))] = city.abb
                split_name = [x + '/' for x in split_name]
                short_name = ''.join(split_name)[0:-1]
        return short_name

    def _deep_search(self, search_value, ignore_case=True, regexp=False):
        
        return t.dict_search(self._data, search_value, ignore_case=ignore_case,
                             regexp=regexp)
                           
    
    def _alias_search(self, query, ignore_case=True, refined=False, 
                      ignore_mutations=False, regexp=True):
        
        '''
        Aliased deep search for the query inside the dictionary. To perform the
        search, it splits the name using / and then deep searches for each 
        element inside the entry. If the search is refined it first tries to 
        guess what is the city name contained in the query and if it is found
        then adds it to the list of things to seach for. 
        
        Then it deep searches for all the split parts of the query inside the entry
        and returns the number of matches. 
        
        if ignore_mutations = False, then the number of matches return 0 if
        the mutations in the query does not exactly match the mutations in the 
        entry.
    
        if regexp = True, comparison is made with regexp which is the defaul 
        behaviour
        
        '''
  
        
        search_names = query.split('/')
        
        # find the mutations in the search query and the entry long name
        # and return 0 matches if they dont match
        if not ignore_mutations:
            mut_pattern = '([ARNDCEQGHILKMFPSTWYV]\d{2,3}[ARNDCEQGHILKMFPSTWYV])'
            long_split = t.flatten_list([z.split('-') for z in t.flatten_list([x.split('_') for x in self.long.split('/')])])
            
            muts1 = t.flatten_list([re.findall(mut_pattern, x) for x in search_names])
            muts2 = t.flatten_list([re.findall(mut_pattern, x) for x in long_split])
            
            if set(muts1) != set(muts2):
                return 0
            
        if refined:   
            
            city_found = False
            found_city_names = []
            for part in search_names:
                if len(city_list.search(part)) > 0: city_found = True
                
            if not city_found:
               
                for part in search_names:
                    
                    found_cities = city_list.search(part)
                    
                    if len(found_cities) > 0:
                        
                        for found_city in found_cities:
                            found_city_names += [found_city.name]
                            city_found = True
                            break
                    
                    elif 'GYRF' in part:
                        found_city_names +=['GYRFALCON', 'WASHINGTON']
                        city_found = True
                        break
                
                if not city_found:
                    
                    for part in search_names:
                    
                        I,dist = city_list.esearch(part)
                        
                        if dist/len(part) <= 0.25:
                            
                            found_city_names += [x.name for x in I]
                            city_found = True
                            break
                

                    
            search_names = t.flatten_list([z.split('-') for z in t.flatten_list([x.split('_') for x in query.split('/')])])
            search_names += found_city_names
        
        matches = 0
         
        for name in search_names:
            
            n_alpha = [x.isalpha() for x in name].count(True)
            n_numeric = [x.isnumeric() for x in name].count(True)
            if len(name)>4 or n_numeric>3 or n_alpha>3:
            
                if self._deep_search(name, ignore_case=ignore_case, regexp=regexp):  
                    matches += 1
        
        
        return matches
     
class SerumEntry(Entry):

        def __init__(self, data):
            super(SerumEntry, self).__init__(data)
            
            if 'strain_id' in data:
                self.strain_id = data['strain_id']
                
            else:
                self.strain_id = None
    
        def __str__(self):
            return 'Serum {} with name {} and strain id {}'.format(self.id, self.long, self.strain_id)
        
        def __repr__(self):
            return '({}, {})'.format(self.id, self.long)
              
class AntigenEntry(Entry):     
    
        def __init__(self, data):
            super(AntigenEntry, self).__init__(data)
            
            if 'parent_id' in data:
                self.parent_id = data['parent_id']
            elif "wildtype" in data and data['wildtype']:
                self.parent_id = self.id 
            else:
                self.parent_id = None
    
        def __str__(self):
            return 'Antigen {} with name {}'.format(self.id, self.long)
        
        def __repr__(self):
            return '({}, {})'.format(self.id, self.long)
              
        
class Dataset:

    '''parent class from which SerumDataset and AntigenDataset classes 
    are constructed. If an id_subset is provided then only those that are in 
    the subset are parsed.'''
    
    def __init__(self, data, id_subset=None): 
        
        assert isinstance(data, list), 'dataset should be a list'
        assert all(isinstance(x, MAPTYPES) for x in data), 'all elements of the dataset should be of map type'
        assert all(['id' in x.keys() for x in data]), 'all elemenys of the dataset should contain the key id'
        
        if id_subset is not None: assert(isinstance(id_subset, list)), 'id_subset should be a list'
        
        self._data = []
        
        if id_subset is None:
            self._data = data
        else:
            for id1 in id_subset:
                I = [x['id'] == id1 for x in data]
                
                if id1 == '':
                    self._data.append({'id':'None', 'long':'None'})
                elif I.count(True) >= 1:
                    self._data += list(compress(data,I))
                    if I.count(True) > 1:
                        print('Warning: more than one data with id {} found inside the dataset'.format(id1))
                elif I.count(True) == 0:
                    print('Warning: no data with id {} found inside the dataset'.format(id1))
                
            
        self.entries = []
        self.ids = []
        self.longs = []

      
    def get_entry(self,  search_value, search_method='id'):
        
        assert search_method in ['id', 'long' ], 'search_method should be id or long'
        search_method = search_method + 's' 
        lookup_list = getattr(self, search_method)
        
        
        I = [index for index, value in enumerate(lookup_list) if len(re.findall(value, search_value))>0 ]
         
        return [self.entries[x] for x in I]
        
       
        
    def deep_search(self, search_value, ignore_case=True):
        
        I = [index for index, entry in enumerate(self.entries) if 
             entry._deep_search(search_value, ignore_case=ignore_case)]
        
        return [self.entries[x] for x in I]
       
    
    def aliased_search(self, search_name, ignore_case=True, regexp=True):
        
        I = np.array([entry._alias_search(
            search_name, ignore_case=ignore_case, 
            regexp=regexp) for entry in self.entries])
        
        max_hits = np.max(I)
        
        if max_hits < 2:  #  if there isn't match with atleast two hits, try a more refined search
            I = np.array([entry._alias_search(
                search_name, ignore_case=ignore_case, 
                refined=True, regexp=regexp) for entry in self.entries])
        
        max_hits = np.max(I)
        
        if max_hits < 2:
            return []
        else:
            return [self.entries[x] for x in np.argwhere(I==max_hits).flatten()]

    def general_health_check(self):
        '''
        1- Datasets are required to have unique ids
        2- Datasets with same long names are required to be non-identical (via deep_equality)
        3- ids are required to be strings
        4- longs are required to be strings

        

        '''

        if len(self.ids) != len(set(self.ids)):
            
            for entry_id in self.ids:
                entries = self.get_entry(entry_id)
                if  len(entries) != 1:
                    
                    print('Entry with id {} was encountered {} times. They are: '.format(entry_id, len(entries)))
                    
                    for nonunique_entry in entries: print(nonunique_entry)
                    
                    print('')
        
        # below we check for entry uniqueness for entries with same long names
        # This is to guard against accidental duplicate entries of the same data.
        
        if len(self.longs) != len(set(self.longs)):
            
            for entry_long in self.longs:
                entries = self.get_entry(entry_long, search_method='long')
                if  len(entries) != 1:
                    
                    for ind1,nonunique_entry1 in enumerate(entries): 
                        for ind2,nonunique_entry2 in enumerate(entries[ind1+1:], start=ind1+1): 
                            if t.deep_eq(nonunique_entry1, nonunique_entry2):
                                
                                print('Following entries with long name {} are identical'.format(entry_long))
                                print(nonunique_entry1)
                                print(nonunique_entry2)
            
                            
        failed_entries = [x for x in self.ids if not isinstance(x, str) ]
        if len(failed_entries)>0:
            print('Following entries have ids which are not strings:')
            for entry in failed_entries: print(*self.get_entry(entry))
            print('')            
        
        failed_entries = [x for x in self.longs if not isinstance(x, str) ]
        if len(failed_entries)>0:
            print('Following entries have longs which are not strings:')
            for entry in failed_entries: print(*self.get_entry(entry, search_method='long'))
            print('')      
        
    def get_all_fields(self):
        
        all_fields = set(['id', 'long'])  # these are guaranteed to exist by
                                             # assertions on construction
        
        for entry in self.entries:
            all_fields = all_fields.union(entry.properties.keys())
        
        return all_fields
    
class SerumDataset(Dataset):

    def __init__(self, data, id_subset=None):
        super(SerumDataset, self).__init__(data, id_subset)
        
        self.entries = [SerumEntry(x) for x in self._data]    
        self.ids = [x.id for x in self.entries]
        self.longs = [x.long for x in self.entries]
        self.strain_ids = [x.strain_id for x in self.entries]
        
        self.general_health_check()
        
    

class AntigenDataset(Dataset):

    def __init__(self, data, id_subset=None):
        super(AntigenDataset, self).__init__(data, id_subset)
        
        self.entries = [AntigenEntry(x) for x in self._data]   
        self.ids = [x.id for x in self.entries]
        self.longs = [x.long for x in self.entries]
        self.parent_ids = [x.parent_id for x in self.entries]
        
        self.general_health_check()


