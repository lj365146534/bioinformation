#########################################
#			ja@2016/9/13				#
#	ref:http://www.afenxi.com/post/5413	#
#			ja@2016/9/13				#
#########################################

library(devtools)
#install_github("badbye/baidumap")
library(baidumap)

# 随便输入几个经纬度坐标
lon = matrix(c(117.93780, 24.55730, 117.93291, 24.57745, 117.23530, 24.64210,
               117.05890, 24.74860), byrow=T, ncol=2)
# 将经纬度坐标转换成真实地理信息
location = getLocation(lon, formatted = T)
location

##   lon=117.9378;lat=24.5573 lon=117.93291;lat=24.57745
## "福建省厦门市海沧区坂南路" "福建省厦门市海沧区大溪路"
##   lon=117.2353;lat=24.6421   lon=117.0589;lat=24.7486
##       "福建省漳州市南靖县"   "福建省漳州市南靖县X607"

# 获取厦门大学经纬度坐标，返回json格式文件
getCoordinate('厦门大学') # json
## 参考文档
##getCoordinate(address, city = NULL, output = "json", formatted = F)
##参数：
##address：地址
##city：可选项，地质所在的城市
##output：json或者xml格式
##formatted：F返回原有的json或者xml格式，而T返回的是经纬度的矩阵
##可以同时多个地点
##getCoordinate(c('北京大学', '清华大学', '人民大学'), formatted = T)  

ad <- getCoordinate('厦门大学', formatted = TRUE)
names(ad) <- NULL

# 绘制地图
# 自己修改了一些参数，并将修改后的package挂在github上，所以我选择从github上安装ggmap包。
# install_github("fibears/ggmap")
library(ggplot2)
library(ggmap)
p <- getBaiduMap("厦门市思明区",zoom = 6,messaging = TRUE)

## Map from URL : http://api.map.baidu.com/staticimage?width=400&height=400¢er=118.13453488213,24.468728076403&zoom=12&scale=2

##		getBaiduMap(location, width = 400, height = 400, zoom = 10, scale = 2,
##      		      color = "color", messaging = TRUE)
##
##		参数：
##		location:包含经度和维度的向量或者是一个矩阵,或者可以是一个字符串表示地址；经纬度和地址将作为地图的中心点
##		width，height：map的宽和高
##		zoom：map的缩放比例，是一个整数，从3（洲）到21（building），默认值是10
##		scale:像素数
##		color："color" or "bw"，表示有色或者是黑白
##		messaging:逻辑语句，决定是否输出下载数据的信息


ggmap(p) + geom_point(aes(x=ad[1], y =ad[2]),color='red',size=5)
 
 

