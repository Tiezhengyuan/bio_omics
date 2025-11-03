import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class PlotSeqRegion:
    def __init__(self, ax, ft:pd.DataFrame=None, chain_len:int=None):
        self.ax = ax
        self.ft = ft
        self.chain_len = chain_len if chain_len else \
            ft[ft['code']=='W']['length'][0]
        self.xlim_max = 100
        self.ylim = 1
        self.height = .3

    def add_dia(self, region_name:str, label:str=None):
        regions = self.ft[self.ft['code']==region_name]
        if len(regions)==0:
            print('wrong region name:', region_name)
            return None

        self.lines= None
        self.boxes = None
        self.baseline = self.ylim
        self.add_lines(regions)
        self.add_boxes(regions)
        label = label if label else region_name
        self.add_label(label)
        self.ylim += 1
        
    def add_lines(self, regions):
        x1 = regions['end'][:-1].tolist()
        x2 = regions['start'][1:].tolist()
        # add line to the head if needed
        if regions['start'].iloc[0] > 0:
            x1.insert(0, 1)
            x2.insert(0, regions['start'].iloc[0])
        # add line to the tail if needed
        if self.chain_len > int(regions['end'].iloc[-1]):
            x1.append(regions['end'].iloc[-1])
            x2.append(self.chain_len)

        # update lines
        self.lines = pd.DataFrame({'x1': x1, 'x2': x2,})
        self.lines['y1'] = self.baseline
        self.lines['y2'] = self.baseline
        for i, row in self.lines.iterrows():
            self.ax.plot(row[:2], row[2:],
                color='black', linewidth=1, linestyle='-')

    def add_boxes(self, regions):
        self.boxes = pd.DataFrame({
            'x1': regions['start'].tolist(),
            'x2': regions['end'].tolist(),
        })
        for i, row in self.boxes.iterrows():
            row = row.tolist()
            # xy, width, height
            rect = patches.Rectangle(
                [row[0], self.baseline - self.height/2],
                row[1] - row[0], self.height,
                linewidth=1, edgecolor='black', facecolor='lightblue'
            )
            self.ax.add_patch(rect)
            self.ax.text(row[0]+2, self.baseline+self.height/2, row[0],
                fontsize=6, ha='right', va='bottom', color='blue')
            self.ax.text(row[1]-2, self.baseline-self.height/2, row[1],
                fontsize=6, ha='left', va='top', color='red')

    def add_label(self, title):
        self.ax.text(-10, self.baseline, title, ha='right', va='center')

    def set_ax(self, x_left:int=60, x_right:int=20):
        xlim_max = self.chain_len + x_right
        self.ax.set_xlim(-x_left, xlim_max)
        self.ax.set_ylim(.6, self.ylim+.5)

        self.ax.axis('off')
        # self.ax.yaxis.set_visible(False)
        # xticks = [0, 10, 100]
        # self.ax.set_xticks(xticks, labels=xticks)