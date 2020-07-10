import os
import pandas as pd
import pastas as ps

import sys
sys.path.insert(1, "../..")

#Werkt niet:
import pastastore as pst

path = "./pystore"
name = "my_second_connector"
conn2 = pst.PystoreConnector(name, path)