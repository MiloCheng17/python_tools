import numpy as np
import pandas as pd
import sys, re

def atoi(text):
        return int(text) if text.isdigit() else text


def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def heatmap(data,row_labels,col_labels,ax=None,cbar_kw={},cbarlabel="",**kwargs):
    if not ax:
        ax = plt.gca()

    # Plot the heatmap    
    im = ax.imshow(data,**kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im,ax=ax,**cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    #cbar.set_clim(-0.5,1)

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",rotation_mode="anchor")

    # Turn spines off and create white grid.
    #ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}", textcolors=("black", "white"),threshold=None, **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2

    # Set default alignment to center, but allow it to be overwritten by textkw.
    kw = dict(horizontalalignment="center",verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


df = pd.read_csv('test/m143-2ts-distmat.csv',header=0)
names = df.columns
#for i in range(len(names)):
#    print(names[i])

num = len(names)
res_unique = np.unique(df['chain-resid'])
res_unique = sorted(res_unique,key=natural_keys)

### Typical hydrogen bond length less than 3.5 A ###

num_res = len(res_unique)
hbond = pd.DataFrame(0,columns=res_unique,index=res_unique)

for i in range(num-2):
    if df['atom'][i][0] == 'H':
        for j in range(2,num):
            if df.iloc[i,j].astype(float) <= 3.5 and df.iloc[i,j] > 0 and df['atom'][j-2][0] != 'H':
                idx1 = res_unique.index(df['chain-resid'][i])
                idx2 = res_unique.index(df['chain-resid'][j-2])
#                print(idx1,idx2)
### The H atom as the donor from residue on the row, H acceptor is the residue on the column
                hbond.iloc[idx1,idx2] += 1
#                print(df.iloc[i,j],df['chain-resid'][i],df['atom'][i],names[j])
                
#print(hbond)                
#
#dfh = df[df['atom'].str.startswith('H')]
#print(dfh.head())
#
#drange = df[(0.0 < df) & (df< 3.5)]
#print(drange)

#dist = df.iloc[:,2:].astype(float)
#print(dist)

#    ### plot distance map ###
fig, ax = plt.subplots()
im, cbar = heatmap(hbond, res_unique, res_unique, ax=ax, cmap ='Spectral',cbarlabel='distance')
texts = annotate_heatmap(im, valfmt="{x:.2f}")

fig.tight_layout()
plt.show()
