'''
'''
import os
import json
import pickle
from datetime import date
import pandas as pd
from .dir import Dir

class Utils:

    @staticmethod
    def count_iterator(iter):
        n = 0
        for _ in iter:
            n += 1
        return n

    @staticmethod
    def update_dict(input:dict, key, value):
        if key not in input:
            input[key] = value
        else:
            if input[key] != value:
                if isinstance(input[key], list):
                    if value not in input[key]:
                        input[key].append(value)
                else:
                    input[key] = [input[key], value]
        return input


    @staticmethod
    def update_nested_dict(input:dict, keys:list, value, current=None):
        if len(keys) == 0:
            return input
        Utils.update_dict(input, keys[0], value)
        Utils.update_nested_dict(input, keys[1:], value, current)

