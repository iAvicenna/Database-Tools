#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 23:16:42 2020

@author: Sina

Test for the package DatabaseTools

"""
import unittest
import sys
import os
import json
sys.path.append('../modules/')
from tools import dict_search, dict_key_search, deep_eq, convert_titer
from CityList import CityList
from TiterTable import TiterTable

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

with HiddenPrints():
    from AntigenicDatabase import AntigenDataset, SerumDataset



class TestTools(unittest.TestCase):
    
    def test_deep_eq(self):
        
        dict1 = {'A':'B', 'C':{'D':{'A':{'A':'DDA1d'}}}}
        dict2 = {'A':'B', 'C':{'D':{'A':{'A':'DDA1d'}}}}
        dict3 = {'A':'B', 'c':{1:{'A':{'A':'DDA1d'},'B':5}}}
        dict4 = {'A':'B', 'c':{1:{'A':{'A':'DDA1d'},'B':5}}}
        dict5 = {}
        dict6 = {}
        
        self.assertTrue(deep_eq(dict1, dict2))
        self.assertFalse(deep_eq(dict1, dict3))
        self.assertTrue(deep_eq(dict3, dict4))
        self.assertTrue(deep_eq(dict5,dict6))
        self.assertFalse(deep_eq(dict1,dict5))
        
        
    def test_dict_key_search1(self):
        
        dict1 = {'A':'B', 'C':{'D':{'a':{'A':'D'}}}}
        
        self.assertTrue(deep_eq(dict_key_search(dict1, 'A', ignore_case=True), [('A','B'), ('a',{'A':'D'}), ('A','D')]))
        
        self.assertTrue(deep_eq(dict_key_search(dict1, 'A'), [('A','B'), ('A','D')]))
        
    def test_dict_key_search2(self):
        
        dict1 = {'A':'B', 'C':{'D':{'a':{'parent_id':'D'}}}}
        
        self.assertTrue(deep_eq(dict_key_search(dict1, 'parent', regexp=True), [('parent_id','D')]))
        self.assertFalse(deep_eq(dict_key_search(dict1, 'parent'), [('parent_id','D')]))

    def test_dict_search(self):
        
        dict1 = {'A':'B', 'C':{'D':{'A':{'A':'DDAd'}}}}
        
        self.assertTrue(dict_search(dict1, 'ddad'))
        self.assertFalse(dict_search(dict1, 'ddad', ignore_case=False))
    
    def test_convert_titer(self):
       
       titers = ['20/40', '<320/640', '640/>1280', '<20', '>5120', '>>5120', 
                 '320.0', '80.0/160', 'A', '<<40/80', '40/A']
       results = []
       expected_results = [28, 320, 1280, 10, 10240, 'ValueError', 320, 113, 
                           'ValueError', 'ValueError', 'ValueError']
       
       for ind,titer in enumerate(titers):
           try:
               with HiddenPrints():
                   results.append(convert_titer(titer))
           except ValueError:
               results.append('ValueError')
        
           self.assertEqual(expected_results[ind], results[ind])
           
class TestCityLists(unittest.TestCase):

    def setUp(self):
        with HiddenPrints():
            self.city_list = CityList()
            
        
    def test_search(self):
        
        query = 'HONGKONG'
        result = self.city_list.search(query)
        self.assertEqual(result[0].name, 'HONGKONG')
    
        query = 'HONGKONG'
        result = self.city_list.search(query, aliasing=True)
        self.assertEqual(set([x.name for x in result]), set(['HONGKONG','HONG-KONG']))
       
        query = 'HONGK'
        result = self.city_list.search(query)
        self.assertEqual(result[0].name, 'HONGKONG')
        
        query = 'HONGK'
        result = self.city_list.search(query, exact_match=True)
        self.assertEqual(len(result), 0)
    
        query = 'HONGKONG/89'
        result = self.city_list.search(query, split=True)
        self.assertEqual(result[0].name, 'HONGKONG')
    
        query = 'HONGKONG/89'
        result = self.city_list.search(query, aliasing=True, split=True)
        self.assertEqual(set([x.name for x in result]), set(['HONGKONG','HONG-KONG']))
    
    
    def test_esearch(self):

        query = 'KINGKONG'
        result = self.city_list.esearch(query)
        self.assertEqual(set([x.name for x in result[0]]), set(['HONGKONG','HONG-KONG']))
        self.assertEqual(result[1],2)
        
        query = 'KINGKONG'
        result = self.city_list.esearch(query, aliasing=False)
        self.assertEqual(set([x.name for x in result[0]]), set(['HONGKONG']))
        self.assertEqual(result[1],2)        
        
    
class TestAntigenDataset(unittest.TestCase):
        
    def setUp(self):
        import os
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        
        with open(script_path + '/test_datasets/test_antigens.json', 'r') as fileobj:
            self.antigens_json = json.load(fileobj)
         
        self.antigen_dataset = AntigenDataset(self.antigens_json)    
        
    def test_antigen_getentry(self):
             
        for entry in self.antigen_dataset.entries:
            entries_found = self.antigen_dataset.get_entry(entry.id)
            self.assertEqual(len(entries_found), 1) 
            self.assertEqual(entries_found[0].long, entry.long) 
        
        
    def test_antigen_deepsearch(self):
        
        entries_found = self.antigen_dataset.deep_search("clade 2.2.x")
        self.assertEqual(len(entries_found), 1) 
        self.assertEqual(entries_found[0].id, "ARTLF7")
        
    def test_antigen_aliasedsearch(self):
        
        entries_found = self.antigen_dataset.aliased_search("A/VIETNAM/1194/2004-NIBRG-14")
        self.assertEqual(len(entries_found), 1)
        self.assertEqual(entries_found[0].id, '14846I')
        
        entries_found = self.antigen_dataset.aliased_search("HONG-KONG/213/2003")
        self.assertEqual(len(entries_found), 1)
        self.assertEqual(entries_found[0].id, "2U7GA8")
        
class TestSerumDataset(unittest.TestCase):
        
    def setUp(self):
        import os
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        
        with open(script_path + '/test_datasets/test_antisera.json', 'r') as fileobj:
            self.sera_json = json.load(fileobj)
         
        self.serum_dataset = SerumDataset(self.sera_json)    
        
    def test_serum_getentry(self):
             
        for entry in self.serum_dataset.entries:
            entries_found = self.serum_dataset.get_entry(entry.id)
            self.assertEqual(len(entries_found), 1) 
            self.assertEqual(entries_found[0].long, entry.long) 
        
    def test_serum_deepsearch(self):
        
        entries_found = self.serum_dataset.deep_search("test sera")
        self.assertEqual(len(entries_found), 1) 
        self.assertEqual(entries_found[0].id, "77POTS")
        
    def test_antigen_aliasedsearch(self):
        
        entries_found = self.serum_dataset.aliased_search("A/BARN-SWALLOW/HONG-KONG/D10-1161/2010_SJ-003")
        self.assertEqual(len(entries_found), 1)
        self.assertEqual(entries_found[0].id, 'CC042E')
        
class TestTiterTable(unittest.TestCase):    
    
    def setUp(self):
        import os
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        
        with open(script_path + '/test_datasets/test_antigens.json', 'r') as fileobj:
            self.antigens_json = json.load(fileobj)
        
        with open(script_path + '/test_datasets/test_antisera.json', 'r') as fileobj:
            self.sera_json = json.load(fileobj)
            
        with open(script_path + '/test_datasets/test_results.json', 'r') as fileobj:
            self.results_json = json.load(fileobj)    
         
           
        
 
    def test_load(self):   
        
        titer_table = TiterTable(self.results_json[0]['results'][0], self.sera_json, self.antigens_json)
 
if __name__ == '__main__':
    suite_tools = unittest.TestLoader().loadTestsFromTestCase(TestTools)
    suite_citylists = unittest.TestLoader().loadTestsFromTestCase(TestCityLists)
    suite_antigens = unittest.TestLoader().loadTestsFromTestCase(TestAntigenDataset)
    suite_sera = unittest.TestLoader().loadTestsFromTestCase(TestSerumDataset)
    suite_titertable = unittest.TestLoader().loadTestsFromTestCase(TestTiterTable)


    unittest.TextTestRunner(verbosity=2).run(suite_titertable)  
    unittest.TextTestRunner(verbosity=2).run(suite_tools)  
    unittest.TextTestRunner(verbosity=2).run(suite_citylists)  
    unittest.TextTestRunner(verbosity=2).run(suite_antigens)  
    unittest.TextTestRunner(verbosity=2).run(suite_sera)       