import os
import os.path
from PyPDF2 import PdfFileReader, PdfFileWriter

import time
def getFileName(filepath):
    list1=[]
    for root, dirs, files in os.walk(filepath,topdown=False):
        for name in files:
            if name.endswith(".pdf"):
                list1.append(os.path.join(root,name))

    return list1

def MergePDF(filepath, outfile):
    output=PdfFileWriter()
    outputPages = 0
    pdf_fileName = getFileName(filepath)
    for each_file in pdf_fileName:
        print("adding %s" % each_file)

        input = PdfFileReader(open(each_file,"rb"))

        if input.isEncrypted == True:
            input.decrypt("map")

        pageCount = input.getNumPages()
        outputPages += pageCount
        print("%s has %d pages" % (each_file, pageCount))

        for iPage in range(pageCount):
            output.addPage(input.getPage(iPage))
        
        output.addBookmark(
            title = each_file.split("\\")[-1],pagenum=outputPages - pageCount)
            
        print("All Pages Number: " + str(outputPages))

        outputStream = open(os.path.join(filepath,outfile),"wb")
        output.write(outputStream)
        outputStream.close()
        print("finished")

if __name__ == '__main__':
    file_dir = ...
    out = u"report.pdf"
    MergePDF(filepath=file_dir,outfile=out)
    