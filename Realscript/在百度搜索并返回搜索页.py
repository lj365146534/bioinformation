import urllib
my_text = '你想搜索的东西'
url = 'http://www.baidu.com/s?wd=%s' % urllib.request.quote(my_text)
result = urllib.request.urlopen(url).read()  # 搜索结果页的内容
print(result)      #可能会看不懂

aa = open('aa.html','wb')    #因为result不是UNICODE，所以要用二进制'b'
aa.write(result)
aa.close()
 
# 想得到网页list,需要解析result,然后从result中把网页提取出来就可以了