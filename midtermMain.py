from pyne import data
from pyne import nucname
import math
import pyne
import warnings
#gets rid of warnings from pyne
from pyne.utils import toggle_warnings
toggle_warnings()

'''Reads peaks from a spectrum and presents likely elements.
Author: Caleb Warren and Michelle Kuchera
'''

def getNum(elementString):
    """Takes in an Element and returns the element number.

    Args:
        elementString: The element in string form.

    Returns:
        the element number associated with this element.

    """    
    num = int(elementString[-1])
    try:
        num = num +int(elementString[-2])*10
    except:
        num = num
    try:
        num += int(elementString[-3])*100
    except:
        num = num
    return num

def countIntensityPeaks(intensities, intensityList,
                        intensityUncertainty):
    """Counts the number of similar intensity peaks.

    Args:
        intensities: all peaks of this element to comare to our data.
        intensityList: The spectrum of intensities that our dataset holds.
        intensityUncertainty: The uncertainties of intensity peaks ( in 
            order) in our dataset
    """    
    #for each true intensity of this element, compare to our table (e)
    counts = 0
    for intense in intensities:
        c = intense[0]
        u = intense[1]
        #catch 'not a number' errors
        if math.isnan(c):
            c = 0
        if math.isnan(u):
            u = 0
        
        if c >0:
            #for each intensity peak in our data...
            for i in range(len(intensityList)):
                #determine if this intensity lies within the intensity range
                if (intensityList[i]-intensityUncertainty[i] <= c+u and 
                    intensityList[i]+intensityUncertainty[i]>=c-u):
                    counts = counts + 1    
    return counts

def findParents(energyList, energyUncertainty):
    """Find parent elements of a given intensity peak.

    Args:
        energyList: The energies with an intensity peak in our dataset.
        energyUncertainty: Energy Uncertainties associated with each
            energy (in order) of our dataset.  
    """    
    elements = []
    for i in range(len(energyList)):
        for j in data.gamma_parent(energyList[i], 
                                   energyUncertainty[i]):
            if j>0:
                try:
                    if  nucname.zzaaam(j) % 2 == 0:
                        elements.append(j)
                except:
                    doNothing = 1
                    raise
    return elements
    
def checkCountsLists(abundance, name, element, bCE, bC,
                     bCN, bCEN,counts):
    """Checks and updates the best lists of counts for elements.

    Args:
        abundance: List of most abundant elements so far.
        name: The name of the element in question.
        element: The element in question.
        bCE: The list of names of the best elements overall, so far.
        bC: List of best counts for all elements.
        bCN: List of best counts for natural elements.
        bCEN: List of names for natural elements for counts.
        counts: number of similar intensity counts for this element.
    """
        #make sure we haven't seen this element before
    isThere = False
    for e in range(len(bCE)-1):
        if bCE[e] == bCE[-1]:
            isThere = True
    smallest = 0
    if not isThere and abundance > 0.0 and abundance < 1.0:
        #update counts lists
        for c in range(len(bC)-1):
            if bC[smallest] > bC[c]:
                smallest = c
        if bC[smallest] < counts:
            bC[smallest] = counts
            bCE[smallest] = bCE[-1]
        if name[-1] != 'M':
            num = getNum(name)
            if num <=90:
                smaller = 0
                for b in range(len(bCN) - 1):
                    if bCN[smaller] > bCN[b]:
                        smaller = b
                if bCN[smaller] < abundance:
                    bCN[smaller] = abundance
                    bCEN[smaller] = bCEN[-1]
                    
def checkAbundanceLists(abundance, name, bA, bAE, bAN, 
                        element, bAEN):
    """Checks and updates the best lists of abundant elements.

    Args:
        abundance: List of most abundant elements so far.
        name: The name of the element in question.
        bA: The list of best overall abundances.
        bAE: List of best abundances names for natural elements.
        element: The element in question.
        bAEN: List of names for natural elements for abundances.
    """
        #make sure we haven't seen this element before
    isThere = False
    for e in range(len(bAE)-1):
        if bAE[e] == bAE[-1]:
            isThere = True
            
    #update abundance lists
    
    if (abundance < 1.0 and not isThere and abundance > 0.0
        and name[-1] != 'M'):
        smaller = 0
        for a in range(len(bA) - 1):
            if bA[smaller] > bA[a]:
                smaller = a
        if bA[smaller] < abundance:
            bA[smaller] = abundance
            bAE[smaller] = bAE[-1]        
            smaller = 0
        #check for natural isotopes
        if data.atomic_mass(element)<=180:
            smaller = 0
            for d in range(len(bAN) - 1):
                if bAN[smaller] > bAN[d]:
                    smaller = d
            if bAN[smaller] < abundance:
                bAN[smaller] = abundance
                bAEN[smaller] = bAEN[-1]  
    
    
def addToLists(bC, bA, bCE, bAE, bCN, bAN, bCEN, bAEN, 
               counts, abundance, name,element):
    """Checks and updates all best lists.

    Args:
        bC: List of best counts for all elements.
        bA: The list of best overall abundances.
        bCE: The list of names of the best elements overall, so far.
        bAE: List of best abundances names for natural elements.
        bCN: List of best counts for natural elements.
        bAN: List of best abundances' names.
        bCEN: List of names for natural elements for counts.
        bAEN: List of names for natural elements for abundances.
        counts: number of similar intensity counts for this element.
        abundance: List of most abundant elements so far.
        name: The name of the element in question.
        element: The element in question.
    """
    if name[-1] != 'M':
        if len(bC) < 11:
            bC.append(counts)
            bA.append(abundance)
            bCE.append(name)
            bAE.append(bCE[-1])
            
        else: 
            bC[-1] = counts
            bA[-1] = abundance
            bCE[-1] = name
            bAE[-1] = bCE[-1]
            
        if len(bCN) < 11 :
            if data.atomic_mass(element)<=180:
                bCN.append(counts)
                bAN.append(abundance)
                bCEN.append(name)
                bAEN.append(bCE[-1])
                
        elif data.atomic_mass(element) <= 180: 
            bCN[-1] = counts
            bAN[-1] = abundance
            bCEN[-1] = name
            bAEN[-1] = bCE[-1]         
    

def presentInfo(bestElement, bAN, bAEN):
    """Presents the best 3 possible elements from this spectrum

    Args:
        bestElement: The best overall element, reguardless of element
            number or abundance.
        bAN: List of best names for counts/abundances.
        bAEN: List of names for natural elements for abundances.
    """
    naturalBestNum= 0.0
    nB = 0    
    num = 0
    bestNum = 0
    print "Possible matches in order of most likely: "
    #go through and get highest abundance, lowest atomic number, 
    #and the best element overall
    bestList = []
    bestNum = getNum(bestElement)
    bestE = ''
    
    for i in range(len(bAEN)):
        if bAN[i] > 0.0:
            b = bAEN[i]
            num = int(b[-1])
            try:
                num = num + int(b[-2])*10
            except:
                num = num
            try:
                num = num+int(b[-3])*100
            except:
                num = num
            if num < bestNum:
                bestNum = num
                bestE = b

    bestList.append(bestE)
    bestList.append(bAEN[nB])
    if bestElement not in bestList:
        bestList.append(bestElement)
    l = ["Natural: ","Heavy/Lab: ","Heavy/Lab: "]
    #print best list
    for e in range(len(bestList)):
        print (l[e] + bestList[e])    
    
def efficiencyAdjust(energyList, energyUncertainty, intensityList,
                     intensityUncertainty):
    """Adjusts efficiency based on Ge detector.

    Args:
       energyList: list of energies in a given spectrum.
       energyUncertainty: Uncertainties for each energy in the spectrum.
       intensityList: List of intensity peaks; one for each energy.
       intensityUncertainty: Uncertainties for each intensity peak in 
       the spectrum.
    """
    for i in range(len(intensityList)):
        #efficient val = 1 for this energy
        y = (-8.45485e-19*energyList[i]*energyList[i]*energyList[i]*
energyList[i]*energyList[i]*energyList[i] + 2.92359e-15*
energyList[i]*energyList[i]*energyList[i]*energyList[i]*energyList[i]
- 3.9719e-12*energyList[i]*energyList[i]*energyList[i]*energyList[i]
+ 2.6714e-9*energyList[i]*energyList[i]*energyList[i] - 9.133059e-7
*energyList[i]*energyList[i] + 0.00014*energyList[i] - 0.00410)
        #HANDLE EFFICIENCY OF UNCERTAINTY LATER!!!
        intensityList[i] = intensityList[i]*intensityList[i]/y
        
def findPeaks(energyList=[], energyUncertainty=[], 
              intensityList=[], intensityUncertainty=[]):
    """Finds relevant 'peak' elements in a spectra.

    Args:
        energyList: energies associated with an intensity peak in our
            dataset.
        energyUncertainty: uncertainties of energy values (in order).
        intensityList: The spectrum of intensities in our dataset. 
        intensityUncertainty: The uncertainties of intensity peaks (in 
            order) of dataset
    """        
    #initialize
    posEl = []
    counts = 0
    bestCount = 0
    bestElement = 'FirstRun'
    abundance = 0.0
    bCE = []
    bAE = []
    bA = []
    bC = []
    bCEN = []
    bAEN = []
    bAN = []
    bCN = []
    
    #for each energy, find the parents associated with that energy
    posEl = findParents(energyList, energyUncertainty)
    
    for element in posEl:
        counts = 0
        #find all intensities associated with this element
        intensities = data.gamma_photon_intensity(element)
        #find unicode name of element
        name = nucname.name(element)
        #count number of intensity peaks that this set has in common
        #with this element
        counts = countIntensityPeaks(intensities, intensityList, 
                                     intensityUncertainty)
                        
        #calculate relative counts and natural abundance     
        counts = 1.0*counts/len(intensities)
        abundance = data.natural_abund(element)
        #update lists
        addToLists(bC, bA, bCE, bAE, bCN, bAN, bCEN, bAEN, 
                   counts, abundance, name,element)
        
        checkAbundanceLists(abundance, name, bA, bAE, bAN, 
                            element, bAEN)                 
        
        checkCountsLists(abundance, name, element, bCE, bC, bCN, 
                         bCEN, counts)
        
        if (counts >bestCount and abundance >0.0 and 
            abundance < 1.0 and name[-1] != 'M'):
            num = getNum(name)
            if bestElement == 'FirstRun':
                bestCount=counts
                bestElement = name
                counts = 0
            
            elif num<getNum(bestElement):
                bestCount=counts
                bestElement = name
                counts = 0
            elif  data.natural_abund(bestElement) < abundance:
                bestCount=counts
                bestElement = name
                counts = 0        
    presentInfo(bestElement, bAN, bAEN)
        
def main():
    """Testing Environment for this program.
    """
    ##Default values for bananas
    #energies in KeV
    #give user option to enter their own data
    userInfo = raw_input("Would you like to enter your own data? (yes or no): ")
    
    energyList = [12.57,16.29,39.43,75.95,92.70,111.95,185.82,
                  238.47,295.20,338.29,351.66,409.87,438.99,510.90,
                  583.19,609.13,727.82,795.25,877.49,903.82,911.15,
                  949.53,968.95,1120.33,1172.87,1459.78,1507.07,
                  1537.02,1586.44,1590.96,1761.70,1843.97]
    
    #energy uncertainty, half of FWHM (KeV)
    energyUncertainty = [1.22,1.22,2.395,1.15,1.11,0.675,1.185,
                         0.81,0.595,0.18,0.535,0.175,0.9,1.46,1.05,1.045,
                         0.14,0.13,0.25,1.285,1.285,1.295,0.73,0.32,0.43,
                         1.345,0.13,0.305,1.35,1.35,1.475,0.455]
    
    #banana intensities, area of each peak in K-areas
    intensityList = [38.1, 27.3, 16.3, 9.68, 3.98, 0.273, 2.58, 2.07, 
                     0.797, 0.339, 1.54, 0.0692, 0.471, 8.27, 1.29, 1.53, 
                     0.463, 0.169, 0.166, 0.205, 1.57, 0.79, 0.641, 0.686, 
                     0.0947, 75.6, 0.0771, 0.015, 0.184, 0.316, 1.04, 0.132]
     
    #uncertainty in area in K-areas
    intensityUncertainty = [0.2018, 0.2053, 0.789, 1.857, 0.406, 
                            0.268, 0.367, 0.329, 0.234, 0.2297, 0.227, 0.137,
                            0.178, 0.241, 0.159, 0.17, 0.149, 0.0917, 0.0885, 
                            0.0497, 0.0594, 0.146, 0.129, 0.154, 0.105, 0.289,
                            0.0376, 0.0299, 0.0224, 0.0248, 0.0577, 0.0339]    
    #adjust for efficiency of detector
    
    findLists = False
    
    #if the user wants to enter their own data, let them
    if userInfo == "yes":
        findLists = True
    
    #while the user hasn't entered a correct format of the data, prompt
    #them to do so
    while findLists:
        
        energyList = raw_input("Please enter the energy of every gamma ray in your spectrum in keV (example: 0.1, 0.5, 100, ...): ")
        energyList = [float(x.strip()) for x in energyList.split(",")]
        
        energyUncertainty = raw_input("Please enter the uncertainty of the energy of every gamma ray in your spectrum in keV: ")
        energyUncertainty = [float(y.strip()) for y in energyUncertainty.split(",")]
        
        intensityList = raw_input("Please enter the intensity of every gamma ray in your spectrum: ")
        intensityList = [float(z.strip()) for z in intensityList.split(",")]
        
        intensityUncertainty = raw_input("Please enter the uncertainty of the intensity of every gamma ray in your spectrum: ")
        intensityUncertainty = [float(w.strip()) for w in intensityUncertainty.split(",")]
        
        firstLen = len(energyUncertainty)
        if (firstLen == len(energyList) and 
            firstLen== len(intensityUncertainty) and 
            firstLen == len(intensityList)):
            print "Processing Data ... "
            findLists = False
        else:
            print "inputs must have the same number of points and be an integer or float only!"
            findLists = True
    efficiencyAdjust(energyList, energyUncertainty, intensityList,
                     intensityUncertainty)
    findPeaks(energyList, energyUncertainty, intensityList,
              intensityUncertainty)
    

if __name__ == "__main__":
    main()
