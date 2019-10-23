file<-'survey_data.csv'

# clean_pp <-function(file){
#   pp<-read.csv(file) #read in data
#   pp<- pp[, c('Site', 'Point', 'Species')] #select relevant columns
#   pp<-pp[-which(pp$Species=='BT'), ] #remove brushtail (BT) data
#   s<- unique(pp$Site) #number of sites
#   p<-unique(pp$Point) #number of points
#   pp_clean<-data.frame(rep(s, each=length(p)), rep(p, times=length(s)), 0) #initialise empty data frame for count data
#   names(pp_clean)<-c('Site', 'Point', 'Count')
#   for (i in s){
#     for (j in p){
#       idx<-which(pp$Site==i & pp$Point ==j)
#       tmp<-pp[idx, ]
#       c<-sum(tmp$Species=='RT') #total number of RT seen at site i, point j
#       #a<-nrow(nrow(tmp[tmp$Species=='None']))
#       if (c!=0){
#       pp_clean[which(pp_clean$Site==i & pp_clean$Point==j), 3]=c
#       } 
#       }
#   }
#   return(pp_clean)
# }

clean_pp<-function(file){
  pp<-read.csv(file) #read in data
  pp<- pp[, c('Site', 'Point', 'Species')] #select relevant columns
  pp$Site_Point <- paste(pp$Site, pp$Point)
  count_RT <- function(species) { sum(species == "RT") }
  counts <- tapply(pp$Species, pp$Site_Point, FUN=count_RT)
  site_points<-names(counts)
  idx<-match(site_points, pp$Site_Point)
  pp_clean<-data.frame(Site=pp$Site[idx], Point=pp$Point[idx], Count=counts)
  pp_clean<-pp_clean[order(pp_clean$Point), ]
  pp_clean<-pp_clean[order(pp_clean$Site), ]
  rownames(pp_clean)<- NULL
}

