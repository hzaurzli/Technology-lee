# import urllib.request
# from bs4 import BeautifulSoup
# import csv
# import pandas as pd
# import time
#
# start = time.time()
# csvFile = '/home/lirz/index.csv'
# data = pd.read_csv(csvFile)
#
# extract = list()
# for geneName in data['gene']:
#     url = 'http://www.uniprot.org/uniprot/?query=' + geneName + '&sort=score'
#     page = urllib.request.urlopen(url)
#     contents = page.read()
#     soup = BeautifulSoup(contents,'html.parser')
#     for tag in soup.find_all('div', class_='content results'):
#         m_body = tag.find('tbody')
#         if m_body is None:
#             print ('no')
#         else:
#             en = m_body.find_all("a")#寻找body中的所有a标签
#             print(en)
#             entry = en[0].text
#
#             organism = en[1].text
#             enn = m_body.find_all("td")#
#             entryName = enn[2].text
#             protein = tag.find('div', class_="protein_names")#精确定位
#             kk = protein.find_all('div', class_="short")
#             proteinName = kk[0].text
#             print(proteinName)
#             extract.append(proteinName + '\n')
#             time.sleep(3)
#
# out = open("/home/lirz/resultFile.txt", 'w')
# out.writelines(extract)
# out.close()
# end = time.time()
# print("hello every body! ")
#
#
#
# print("The result number {}, time cost {}".format(len(extract), end - start))
# '''
#             with open('/home/lirz/resultFile', 'a', newline='') as out:
#                 csv_write = csv.writer(out, dialect='excel')
#                 print(out)
#                 # csv_write.writerow(label)
#                 # csv_write.writerow(result)
# '''
# from bs4 import BeautifulSoup
# import requests
# import re
# import time
# print("Begin!")
# def GetMyGO(url):
#     opline=[]
#     try:
#         mycontent=requests.get("https://www.uniprot.org/uniprot/"+url,
#                    headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) '
#                          'Chrome/67.0.3396.79 Safari/537.36'}).content
#         soup=BeautifulSoup(mycontent,'html.parser', from_encoding='utf-8')
#         GOs=soup.find_all(class_=re.compile("(.+?cellular_component)|(.+?molecular_function)|(.+?biological_process)"))
#         if GOs:
#             for element in GOs:
#                 function="\t".join(element['class'])
#                 for iterm in element.find_all('a'):
#                     GO_number=iterm["href"].split("/")[-1]
#                     name=iterm.text
#                     oplines="**".join([url,function,name,GO_number])+"\n"
#                     opline.append(oplines)
#                     print(oplines)
#         else:
#             opline.append("**".join([url,"NotFound"])+"\n")
#     except:
#         opline.append("**".join([url,"Unknown Wrong"])+"\n")
#     return opline
# workpath="/share/users_home/dgwei/lizhen/plantvirus_20190320/SequenceAlignment/StrainsAlign/"
# opfile=open(workpath+"scrach_GO",'w')
# with open(workpath+"uniprot_list.txt",'r') as F:
#     for line in F:
#         opfile.write(GetMyGO(line.rstrip()))
#         opfile.write("\n")
#         time.sleep(2)
# opfile.close()



import urllib.request
from bs4 import BeautifulSoup
import csv
import pandas as pd
import time

start = time.time()
csvFile = '/home/lirz/index.csv'
data = pd.read_csv(csvFile)

#extract = list()
for geneName in data['gene']:
    url = 'http://www.uniprot.org/uniprot/?query=' + geneName + '&sort=score'
    page = urllib.request.urlopen(url)
    contents = page.read()
    soup = BeautifulSoup(contents, 'html.parser')
    for tag in soup.find_all('div', class_='content results'):
        m_body = tag.find('tbody')
        if m_body is None: #判断有没有
            print('no')
        else:
            en = m_body.find_all("a")  # 寻找body中的所有a标签
            entry = en[0].text
            url1 = 'https://www.uniprot.org/uniprot/'+entry

            page1 = urllib.request.urlopen(url1)
            contents1 = page1.read()
            soup1 = BeautifulSoup(contents1, 'html.parser')

            for tag1 in soup1.find_all('ul', class_='noNumbering molecular_function'):  #区分不同功能MF，CC，BP
                print(tag1)
                en1 = tag1.find_all("a")
                print(en1)
                enn1 = en1[0].text
                print(enn1)
                extract.append(geneName+ enn1 + '\n')
                time.sleep(3)

out = open("/home/lirz/resultFile.txt", 'w')
out.writelines(extract)
out.close()
# end = time.time()