import re
import numpy as np
import collections

MAPTYPES = (dict, collections.abc.Mapping)

    
def dict_search(dictionary, search_value, ignore_case=True, regexp=False):
    '''
    This is a deep search of a value (of str type) inside a dictionary. It also 
    recursively searches inside values that are dictionaries/lists. The search 
    supports regular expressions. 
    
    Parameters
    ----------
    dictionary : dictionary or mapping
        
    search_value : string
    
    ignore_case: If the search_value is a string, whether to ignore case or not
    
    if regexp = True, comparison is made assuming a regular expression
   
    Returns
    -------
    bool
        
    '''
    
    assert isinstance(search_value, str)
    
    if regexp:
        if ignore_case:
            search_expression = re.compile(search_value, re.IGNORECASE)
        else:
            search_expression = re.compile(search_value)
            
        def comparitor(x,y): return x.search(y)
        
    else:
        
        if ignore_case:
            search_expression = search_value.upper()
        else:
            search_expression = search_value
    
        def comparitor(x,y): return x==y
  
    
    for key in dictionary:
        lookup_values = dictionary[key]
        
        if isinstance(lookup_values, list):
            pass
        else:
            lookup_values = [lookup_values]
        
        
        for lookup_value in lookup_values:
            if ignore_case:
                
                try:
                    lookup_expression = lookup_value.upper()
                except AttributeError:
                    lookup_expression = lookup_value
            else:
                
                lookup_expression = lookup_value
                    
            if isinstance(lookup_value, MAPTYPES):  
                if dict_search(lookup_expression, search_value, 
                               ignore_case=ignore_case, regexp=regexp): return True
            
            if isinstance(lookup_value, str):
                if comparitor(search_expression, lookup_expression): return True
                
            
        
    return False


def dict_key_search(dictionary, search_key, ignore_case = False, regexp=False):
    '''
    This is a deep search of a key (of string type) inside a dictionary. 
    It also recursively searches inside values that are dictionaries. The search 
    supports regular expressions and returns a found key value pair. If there 
    are multiple such keys, it returns all of them.
    
    Parameters
    ----------
    dictionary : dictionary or mapping
        
    search_key : string
    
    ignore_case: If the search_value is a string, whether to ignore case or not
   
    Returns
    -------
    list
        
    '''
    
    assert isinstance(search_key, str)
        
    if regexp:
        if ignore_case:
            search_expression = re.compile(search_key, re.IGNORECASE)
        else:
            search_expression = re.compile(search_key)
            
        def comparitor(x,y): return x.search(y)
        
    else:
        
        if ignore_case:
            search_expression = search_key.upper()
        else:
            search_expression = search_key
    
        def comparitor(x,y): return x==y
  
  
    found_keys = []
  
    for key in dictionary:
        lookup_value = dictionary[key]
        
        if ignore_case:
            try:
                key_expression = key.upper()
            except AttributeError:
                key_expression = key
        else:
            key_expression = key
           
        if isinstance(key, str):
            if comparitor(search_expression, key_expression): found_keys += [(key, lookup_value)]
        
        if isinstance(lookup_value, MAPTYPES):  
            found_keys += flatten_list([dict_key_search(
                lookup_value, search_key, ignore_case=ignore_case,
                regexp=regexp)])
        
        
    return found_keys

def deep_eq(_v1, _v2):
  """
  Tests for deep equality between two python data structures recursing 
  into sub-structures if necessary. Works with all python types including
  iterators and generators. This function was dreampt up to test API responses
  but could be used for anything. Be careful. With deeply nested structures
  you may blow the stack.
  
  Information: from samuraisam/deep_eq.py github, modified for compatibility
  to python 3
  """
  import operator
  
  def _deep_dict_eq(d1, d2):
    k1 = sorted(d1.keys())
    k2 = sorted(d2.keys())
    if k1 != k2: # keys should be exactly equal
      return False
    return sum(deep_eq(d1[k], d2[k]) for k in k1) == len(k1)
  
  def _deep_iter_eq(l1, l2):
    if len(l1) != len(l2):
      return False
    return sum(deep_eq(v1, v2) for v1, v2 in zip(l1, l2)) == len(l1)
  
  op = operator.eq
  c1, c2 = (_v1, _v2)
  
  # guard against strings because they are also iterable
  # and will consistently cause a RuntimeError (maximum recursion limit reached)
  for t in [str]:
    if isinstance(_v1, t):
      break
  else:
    if isinstance(_v1, dict):
      op = _deep_dict_eq
    else:
      try:
        c1, c2 = (list(iter(_v1)), list(iter(_v2)))
      except TypeError:
        c1, c2 = _v1, _v2
      else:
        op = _deep_iter_eq
  
  return op(c1, c2)       

def flatten_list(list_to_flatten):
   
    assert all(isinstance(x, list) for x in list_to_flatten), 'All elements of list to be flattened should be lists'
    return [item for sublist in list_to_flatten for item in sublist]

def convert_titer(titer):
    
    '''
    If titer is of the form t1/t2 return geometric mean of t1, t2
    If titer contants <, halve it, if contains > double it (also applies
    to t1, t2 above by recursion)
    
    If titer is not a string, return it as it is. 
    
    Do not allow 0 titers.
    
    Returns the titer as a numerical value, not a string.
    
    '''
    
    
    if isinstance(titer,str):
        
        
        if '/' in titer:    
            vals = titer.split('/')
            
            val1 = convert_titer(vals[0])
            val2 = convert_titer(vals[1])
            
            assert val1 != 0 and val2 != 0, 'Zero titer detected in entry {}'.format(titer)  
            return int(np.sqrt(val1 * val2)) # geometric mean of titers = arithmetic mean of log titers
        
        if '<' in titer:
            try:
                float(titer[1:])
            except ValueError:
                print('Titer {} is not convertible to numerical value.'.format(titer))
                
            val = int(float(titer[1:]) / 2) 
            assert val != 0, 'Zero titer detected in entry {}'.format(titer)
            return val
        
        if '>' in titer:
            try:
                float(titer[1:])
            except ValueError:
                print('Titer {} is not convertible to numerical value.'.format(titer))
                    
            val = int(float(titer[1:]) * 2)  
            assert val != 0, 'Zero titer detected in entry {}'.format(titer)
            return val
        
        try:
            float(titer)
        except ValueError:
            print('Titer {} is not convertible to numerical value.'.format(titer))
                
        return int(float(titer))

    else:
        print('Warning titer is not of a string type, returning it as it is.')
        return titer
        

