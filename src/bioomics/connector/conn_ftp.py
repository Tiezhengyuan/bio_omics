"""
connect FTP using object methods

Single FTP connector of bioomics, merged from the former conn_ftp.py (ConnFTP)
and conn_ftplib.py (ConnFTPlib):
1. object methods keep the download context: url, overwrite, gunzip
2. every file is downloaded through its own short-lived connection
   avoiding timeouts on long-lived logins
3. connection and download failures are retried
   (default 3 times with a 2s backoff between attempts)
4. wrong FTP paths (error_perm 550) are reported without retrying
"""
from pathlib import Path
from ftplib import FTP, error_perm, all_errors
import os
import re
import time
from sh import gunzip


class ConnFTP:

    def __init__(self, url:str, overwrite:bool=None, run_gunzip:bool=None,
            try_times:int=3):
        '''
        arg url: FTP host, anonymous login
        arg overwrite: download files again even when they exist locally,
            default False: existing local files are skipped
        arg run_gunzip: decompress downloaded .gz files
        arg try_times: times a failed connection or download is retried
        '''
        self.url = url
        self.overwrite = bool(overwrite)
        self.run_gunzip = bool(run_gunzip)
        self.try_times = try_times if try_times else 3

    # connection ------------------------------------------------------- #

    def connect(self, endpoint:str=None):
        '''
        connect, login and change directory of the FTP host
        transient failures are retried self.try_times times (2s backoff),
        a wrong FTP path raises immediately without retry
        return ftplib.FTP
        raise IOError when the connection finally fails
        '''
        error = None
        for i in range(1, self.try_times + 1):
            ftp = None
            try:
                ftp = FTP(self.url, timeout=100)
                ftp.login()
                if endpoint:
                    ftp.cwd(endpoint)
                return ftp
            except error_perm as e:
                # wrong FTP path: retrying wouldn't help
                self._quit(ftp)
                raise IOError(f"Failure: FTP path is wrong, "
                    f"{self.url}/{endpoint}. error={e}") from e
            except all_errors as e:
                error = e
                self._quit(ftp)
                if i < self.try_times:
                    time.sleep(2)
        raise IOError(f"Failure: can't connect {self.url}/{endpoint}. " +
            f"Try {self.try_times} times. error={error}")

    @staticmethod
    def _quit(ftp):
        '''close the connection, ignoring errors of a broken connection'''
        if ftp is None:
            return
        try:
            ftp.quit()
        except Exception:
            pass

    @staticmethod
    def _nlst(ftp):
        '''
        list names of the current FTP directory
        return [] for an empty directory (nlst raises error_perm there)
        '''
        try:
            return ftp.nlst()
        except error_perm:
            return []

    def _is_dir(self, ftp, name:str):
        try:
            origin_dir = ftp.pwd()
            ftp.cwd(name)
            ftp.cwd(origin_dir)
            return True
        except Exception as e:
            pass
        return False

    # listing ---------------------------------------------------------- #

    def list_files(self, endpoint:str=None, match:str=None):
        '''
        list ftp path of files
        not recursive
        return a list of tuples (endpoint, file_name)
        '''
        # connect FTP
        ftp = self.connect(endpoint)
        try:
            # scan
            ftp_files = []
            for file_name in self._nlst(ftp):
                if not self._is_dir(ftp, file_name):
                    if match is None or re.search(match, file_name):
                        ftp_files.append((endpoint, file_name))
            return ftp_files
        finally:
            self._quit(ftp)

    # download --------------------------------------------------------- #

    def download_file(self, endpoint:str, file_name:str, local_path:str):
        '''
        download one file from FTP
        one download one connection avoiding timeout,
        connection failures are retried self.try_times times with 2s backoff,
        a wrong FTP path or file name is not retried
        return the local file path (decompressed when run_gunzip),
        None on failure
        '''
        Path(local_path).mkdir(parents=True, exist_ok=True)
        local_file = os.path.join(local_path, file_name)
        unzip_file = local_file.replace('.gz', '')
        # doesn't download if file exists and overwrite is False
        if self.overwrite is False:
            if os.path.isfile(unzip_file):
                return unzip_file
            elif os.path.isfile(local_file):
                return local_file

        # download: a whole connect + retrieve cycle per attempt
        error = None
        for i in range(1, self.try_times + 1):
            ftp = None
            try:
                ftp = FTP(self.url, timeout=100)
                ftp.login()
                if endpoint:
                    ftp.cwd(endpoint)
                with open(local_file, 'wb') as f:
                    ftp.retrbinary(f"RETR {file_name}", f.write)
                print(f"Download {self.url}/{endpoint or ''}/{file_name}")
                error = None
                break
            except error_perm as e:
                # wrong path or file name: retrying wouldn't help
                error = e
                print(f"Failure: can't download data from FTP, the path "
                    f"{self.url}/{endpoint or ''}/{file_name} is wrong. "
                    f"error={e}")
                break
            except all_errors as e:
                error = e
                if i < self.try_times:
                    time.sleep(2)
                else:
                    print(f"Failure: can't download data from FTP, "
                        f"{self.url}/{endpoint or ''}/{file_name}. " +
                        f"Try {self.try_times} times. error={e}")
            finally:
                self._quit(ftp)

        # remove the partial file of a failed download
        if error is not None:
            if os.path.isfile(local_file):
                os.remove(local_file)
            return None

        # unzip .gz file
        if self.run_gunzip and local_file.endswith('.gz'):
            print(f"decompress {local_file} to {unzip_file}")
            gunzip(local_file, '-f')
            return unzip_file
        return local_file

    def download_files(self, endpoint:str=None, match:str=None, local_path:str=None):
        '''
        download files from FTP path
        That is not recursive
        '''
        # list all files
        ftp_files = self.list_files(endpoint, match)
        print(endpoint, ftp_files)

        # download files
        local_files = []
        for current_endpoint, file_name in ftp_files:
            local_file = self.download_file(
                endpoint = current_endpoint,
                file_name = file_name,
                local_path = local_path,
            )
            if local_file and local_file not in local_files:
                local_files.append(local_file)
        return local_files

    def download_tree(self, local_path:str, endpoint:str=None, match:str=None):
        '''
        arg: local_name is determined by os.path.join()
        Download FTP directory recursively
        The files of the endpoint itself are stored in local_path,
        sub-directories are mirrored under local_path
        Each directory is scanned through its own short-lived connection
        and every file is downloaded on its own connection with retries
        '''
        # initialize endpoint
        if endpoint:
            endpoint = endpoint.rstrip('/')
        local_files = []
        scanned = set()

        # scan ftp path: one short-lived connection per directory
        pool = [(endpoint, local_path)]
        while pool:
            current_endpoint, _local_path = pool.pop(0)
            if current_endpoint in scanned:
                continue
            scanned.add(current_endpoint)
            Path(_local_path).mkdir(parents=True, exist_ok=True)

            dirs, files = [], []
            try:
                ftp = self.connect(current_endpoint)
                try:
                    for name in self._nlst(ftp):
                        if self._is_dir(ftp, name):
                            dirs.append(name)
                        elif match is None or re.search(match, name):
                            files.append(name)
                finally:
                    self._quit(ftp)
            except IOError as e:
                print(f"Failure: skip FTP path. {e}")
                continue
            print(current_endpoint, files)

            # queue sub-directories
            for name in dirs:
                sub_endpoint = f"{current_endpoint}/{name}" \
                    if current_endpoint else name
                pool.append((sub_endpoint, os.path.join(_local_path, name)))

            # download files of the current directory
            for file_name in files:
                local_file = self.download_file(
                    endpoint = current_endpoint,
                    file_name = file_name,
                    local_path = _local_path,
                )
                if local_file and local_file not in local_files:
                    local_files.append(local_file)
        return local_files