from __future__ import division
from collections import deque
from decimal import Decimal
from math import log
import sys
import itertools

class Base:
    def __init__(self, nucl=None):
       self.nucl = nucl
       self.next = None

    def __str__(self):
        """Compute the string representation of the base"""
        return "%s" % (self.nucl)

    def __repr__(self):
        """Compute the string representation of the base"""
        return "%s(%s)" % (self.__class__.__name__, self.nucl)
 
class Codon:
    def __init__(self, nucl=None):
        self.head = Base(nucl)
        self.tail = self.head
        self.push(nucl)
        self.push(nucl)

    def add(self, nucl):
        self.push(nucl)
        return self.pull()
 
    def push(self,  nucl=None):
        self.tail.next = Base(nucl)
        self.tail = self.tail.next

    def pull(self):
        pulled_base = self.head
        self.head = self.head.next
        return pulled_base
 
    def pop(self):
        #popped = self.tail
        #self.tail = self.tail.next
        #return popped
        pass

    def __str__(self):
        return "%s%s%s" % (self.head, self.head.next, self.tail)

    def __repr__(self):
        """Compute the string representation of the base"""
        return "%s(%s %s %s)" % (self.__class__.__name__,
                                 repr(self.head),
                                 repr(self.head.next),
                                 repr(self.tail)
                                 )

