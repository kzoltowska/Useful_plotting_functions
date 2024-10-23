import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

def volcano_plot(df, logFC, adjpval, pvalcut, FCcut, label, label_num=10):
    """
    Function to generate volcano plot from differential expression data
    param: df - dataframe with differential expression analysis output
    param: logFC - column name with the logFC values
    param: adjpval - column name with the adjpval
    param: pvalcut - threshold for adj.P.Val
    param: FCcut - threshold for logFC
    param: label - column name with gene names
    param: label_num - [int, "All"] how many labels to plot, default 10 top logFC.

    """
    
    fig=plt.figure(dpi=600)
    ax = plt.axes()
    ax.set_facecolor("white")
    plt.scatter(x=df[logFC],y=df[adjpval].apply(lambda x:-np.log10(x)),s=3,label="Not significant", color=sns.color_palette("dark")[8], alpha=0.6)
    # highlight down- or up- regulated genes
    down = df[(df[logFC]<=-FCcut)&(df[adjpval]<=pvalcut)]
    up = df[(df[logFC]>=FCcut)&(df[adjpval]<=pvalcut)]

    plt.scatter(x=down[logFC],y=down[adjpval].apply(lambda x:-np.log10(x)),s=5,label="Down-regulated",color=sns.color_palette("dark")[3], alpha=0.8)
    plt.scatter(x=up[logFC],y=up[adjpval].apply(lambda x:-np.log10(x)),s=5,label="Up-regulated",color=sns.color_palette("dark")[4], alpha=0.8)

    
    if label_num=="All_sign":
        for index, row in df.iterrows():
            if (abs(row[logFC])>1) and row[adjpval]<0.05:
                plt.text(x=row[logFC],y=-np.log10(row[adjpval]),s=row[label], fontfamily="sans-serif", fontsize=5)
    elif str(label_num).isnumeric():
        df["logFCabs"]=abs(df[logFC])
        df1=df.sort_values(by="logFCabs", ascending=False).copy()
        df1=df1[df1[adjpval]<0.05]
        df1=df1.reset_index()
        for row in range(label_num):
            plt.text(x=df1[logFC][row], y=-np.log10(df1[adjpval][row]), s=df1[label][row], fontfamily="sans-serif", fontsize=5)
    
    plt.xlabel("logFC", size=10)
    plt.ylabel("-logFDR", size=10)
    plt.axvline(-FCcut,color="grey",linestyle="--")
    plt.axvline(FCcut,color="grey",linestyle="--")
    plt.axhline(-np.log10(pvalcut),color="grey",linestyle="--")
    plt.legend()
    plt.show()
