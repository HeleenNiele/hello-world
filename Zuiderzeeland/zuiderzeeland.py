import numpy as np
import matplotlib.pyplot as plt
import easygui
import tkinter
import csv
import pandas as pd

if __name__ == "__main__":
    # print("Geef postcodetabel op.")
    # input_file = easygui.fileopenbox(title="Selecteer de postcodetabel:")

    # print("Geef uitvoermap op.")
    # root = tkinter.Tk()
    # root.withdraw()
    # output_path = tkinter.filedialog.askdirectory(parent=root,initialdir="/",title='Geef een uitvoermap op:')

    postcode_path = r"C:\Users\hnx\Desktop\Projecten\2020\Zuiderzeeland\postcodetabel.xlsx"
    data_path = r"C:\Users\hnx\Desktop\Projecten\2020\Zuiderzeeland\Overzicht bedrijfsafvalwaterlozingen en belasting rwzi 2016 2017 2018.xlsx"
    output_path = r"C:\Users\hnx\Desktop\Projecten\2020\Zuiderzeeland"

    print("postcodes lezen ...")
    postcodes = pd.read_excel(postcode_path,index_col='PostcodeID',usecols=['PostcodeID','Plaats','Gemeente','Latitude','Longitude'])
    print(postcodes.columns)
    