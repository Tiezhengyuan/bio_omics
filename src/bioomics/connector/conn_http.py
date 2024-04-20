'''
download data using HTTP
'''
import os
from requests.models import Response
import urllib3

class ConnHTTP:
    def __init__(self, url:str, local_path:str, overwrite:bool=None):
        self.url = url
        self.local_path = local_path
        self.overwrite = True if overwrite else False
    
    def download_file(self, end_point:str, file_name:str=None):
        if file_name is None:
            file_name = end_point
        try:
            url_path = f"{self.url}{end_point}"
            local_file = os.path.join(self.local_path, file_name)
            res = self.to_file(url_path, local_file)
            return res, local_file
        except Exception as e:
            print(e)
            res = Response()
            res.status = 501
            res.error =e 
            return res, None
    
    def to_file(self, url_path, local_file):
        if self.overwrite or (not os.path.isfile(local_file)):
            res = urllib3.request('GET', url_path)
            if res.status in (200, '200'):
                with open(local_file, 'wb') as f:
                    f.write(res.data)
            return res