import os
import re
import fnmatch

def atoi(text):
    return int(text) if text.isdigit() else text

def nat_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def natSort(unsorted):
    unsorted.sort(key=nat_keys)

def applyFilter(unfiltered, filter):
    return fnmatch.filter(unfiltered, filter)

class Folder:
    def __init__(self, path):
        self.path = path
        self.dirs = []
        self.names = []
        self.files = []
        self.__scoutPath()

    def __scoutPath(self):
        _, self.dirs, self.names = os.walk(self.path).__next__()
        natSort(self.dirs)
        natSort(self.names)
        self.files = [self.path + f for f in self.names]

    def filterDirs(self, filter):
        self.dirs = applyFilter(self.dirs, filter)

    def filterFiles(self, filter):
        self.names = applyFilter(self.names, filter)
        self.files = [self.path + f for f in self.names]
