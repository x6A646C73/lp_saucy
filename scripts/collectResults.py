import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def show_boxplot( data, labels, codeVal, top ):
    """
    Thanks Josh Hemann for the example
    """

    numDists = len( labels )

    fig, ax1 = plt.subplots(figsize=(10,6))
    fig.canvas.set_window_title('Code times ' + str(codeVal))
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5, showfliers=False)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    #ax1.set_title('Comparison of IID Bootstrap Resampling Across Five Distributions')
    ax1.set_xlabel('Integer Program')
    ax1.set_ylabel('Time (ms)')

    # Now fill the boxes with desired colors
    boxColors = ['darkkhaki','royalblue']
    numBoxes = numDists*2
    medians = range(numBoxes)
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX,boxY)
        # Alternate between Dark Khaki and Royal Blue
        k = i % 2
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(data[i])],
                   color='w', marker='*', markeredgecolor='k')

    # Set the axes ranges and axes labels
    ax1.set_xlim(0.5, numBoxes+0.5)
    bottom, top = ax1.get_ylim()
    top += 1
    #bottom += -1
    ax1.set_ylim(bottom, top)
    xtickNames = plt.setp(ax1, xticklabels=np.repeat(labels, 2))
    plt.setp(xtickNames, rotation=45, fontsize=8)

    # Due to the Y-axis scale being different across samples, it can be
    # hard to compare differences in medians across the samples. Add upper
    # X-axis tick labels with the sample medians to aid in comparison
    # (just use two decimal places of precision)
    pos = np.arange(numBoxes)+1
    upperLabels = [str(np.round(s, 2)) for s in medians]
    weights = ['bold', 'semibold']
    for tick,label in zip(range(numBoxes),ax1.get_xticklabels()):
        k = tick % 2
        ax1.text(pos[tick], top-(top*0.05), upperLabels[tick],
                horizontalalignment='center', size='x-small', weight=weights[k],
                color=boxColors[k])

    # Finally, add a basic legend
    plt.figtext(0.80, 0.08,  '     Original Saucy' ,
                backgroundcolor=boxColors[0], color='black', weight='roman', size='x-small')
    plt.figtext(0.80, 0.045, '     Modified Saucy',
                backgroundcolor=boxColors[1], color='black', weight='roman', size='x-small')
    plt.figtext(0.80, 0.015, '*', color='white', backgroundcolor='silver', weight='roman', size='medium')
    plt.figtext(0.815, 0.013, ' Average Value', color='black', weight='roman', size='x-small')

    plt.show()


def main():
    files = os.listdir( "../graph" )

    fCount = 0
    #data = [[],[],[],[],[],[],[],[]]
    data = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    labels = []
    skip = [["dano3mip_c","mitre_c","mkc_c","rentacar_c","arki001_c","p2756_c","swath_c","mod011_c"]]
    #skip = [["dano3mip_c","mitre_c","mkc_c","rentacar_c"], ["arki001_c","p2756_c","swath_c","mod011_c"]]
    
    codeVal = 0

    maxV = 0
    for f in files:
        if f == "cap6000_c":
            continue
        if f == "simple":
            continue
        #if f in skip[0] or f in skip[1]:
        if f in skip[0]:
            continue

        print "Collecting results for %s..." % f,

        saucyFile = open( "../res/"+f+".saucy", "r" )
        saucy2_5File = open( "../res/"+f+".saucy2-5", "r" )

        count = 0
        for line in saucyFile:
            junk = 10e3 * float(line)
            if junk > maxV:
                maxV = junk
            data[2*fCount].append( junk )
            count += 1
            if count >= 1000:
                break

        count = 0
        for line in saucy2_5File:
            if count < 1000:
                junk = 10e3 * float(line)
                if junk > maxV:
                    maxV = junk
                data[2*fCount + 1].append( junk )
                count += 1
            else:
                temp = line.split()
                if temp[0] == "different":
                    weights = temp[3]
                    labels.append( f + " (" + weights + ")" )
                    break

        saucyFile.close()
        saucy2_5File.close()

        fCount += 1
        if fCount == 8:
            codeVal += 1
            show_boxplot( data, labels, codeVal, maxV )
            fCount = 0
            maxV = 0
            #data = [[],[],[],[],[],[],[],[]]
            data = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            labels = []

        print weights, "... Done"

    #data = [[],[],[],[],[],[],[],[]]
    data = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    labels = []
    fCount = 0
    maxV = 0
    for f in skip[0]:
        print "Collecting results for %s..." % f,

        saucyFile = open( "../res/"+f+".saucy", "r" )
        saucy2_5File = open( "../res/"+f+".saucy2-5", "r" )

        count = 0
        for line in saucyFile:
            junk = 10e3 * float(line)
            if junk > maxV:
                maxV = junk
            data[2*fCount].append( junk )
            count += 1
            if count >= 1000:
                break

        count = 0
        for line in saucy2_5File:
            if count < 1000:
                junk = 10e3 * float(line)
                if junk > maxV:
                    maxV = junk
                data[2*fCount + 1].append( junk )
                count += 1
            else:
                temp = line.split()
                if temp[0] == "different":
                    weights = temp[3]
                    labels.append( f + " (" + weights + ")" )
                    break

        saucyFile.close()
        saucy2_5File.close()

        fCount += 1
        if fCount == 8:
            codeVal += 1
            show_boxplot( data, labels, codeVal, maxV )
            fCount = 0
            maxV = 0
            #data = [[],[],[],[],[],[],[],[]]
            data = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            labels = []

        print weights, "... Done"

    '''
    data = [[],[],[],[],[],[],[],[]]
    labels = []
    fCount = 0
    maxV = 0
    for f in skip[1]:
        print "Collecting results for %s..." % f,

        saucyFile = open( "../res/"+f+".saucy", "r" )
        saucy2_5File = open( "../res/"+f+".saucy2-5", "r" )

        count = 0
        for line in saucyFile:
            junk = 10e3 * float(line)
            if junk > maxV:
                maxV = junk
            data[2*fCount].append( junk )
            count += 1
            if count >= 1000:
                break

        count = 0
        for line in saucy2_5File:
            if count < 1000:
                junk = 10e3 * float(line)
                if junk > maxV:
                    maxV = junk
                data[2*fCount + 1].append( junk )
                count += 1
            else:
                temp = line.split()
                if temp[0] == "different":
                    weights = temp[3]
                    labels.append( f + " (" + weights + ")" )
                    break

        saucyFile.close()
        saucy2_5File.close()

        fCount += 1
        if fCount == 4:
            codeVal += 1
            show_boxplot( data, labels, codeVal, maxV )
            fCount = 0
            maxV = 0
            data = [[],[],[],[],[],[],[],[]]
            labels = []

        print weights, "... Done"
    '''


if __name__ == "__main__":
    main()
