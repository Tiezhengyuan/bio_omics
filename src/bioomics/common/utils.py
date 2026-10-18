'''
'''
import os
import socket
import subprocess

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
    def update_nested_dict(input:dict, keys:list, value):
        '''
        set value at the keys path, creating intermediate dicts;
        the leaf merge follows Utils.update_dict semantics.
        NB: the old version wrote every key flat into input and
        returned None
        '''
        if not keys:
            return input
        d = input
        for key in keys[:-1]:
            nxt = d.get(key)
            if not isinstance(nxt, dict):
                nxt = {}
                d[key] = nxt
            d = nxt
        Utils.update_dict(d, keys[-1], value)
        return input

    @staticmethod
    def copy_data(outdir) -> bool:
        data_hostname = os.getenv('DATA_HOSTNAME')
        if not data_hostname:
            print("No data hostname is assigned, skip synchronization.")
            return False
        
        execute_hostname = socket.gethostname()
        if data_hostname == execute_hostname:
            print(f"Skip synchronization, identical IP={data_hostname}")
        else:
            print(f"Try to copy {outdir} from {execute_hostname} to {data_hostname}.")
            data_host_user = os.getenv('DATA_HOST_USER')
            data_host_ip = os.getenv('DATA_HOST_IP')
            cmd = [
                'rsync', '-avz', '--progress',
                '--mkpath', '--remove-source-files',
                f"{outdir}/*",
                f"{data_host_user}@{data_hostname}:{outdir}/",
            ]
            try:
                cmd = ' '.join(cmd)
                print(cmd)
                res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
                print('Synchronize results: ', res.stdout, res.stderr)
                return True
            except Exception as e:
                print(f"ERROR: Synchronization between {outdir}, err={e}")
        return False