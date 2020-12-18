#  Animate multiple variables of hourly ERA5-Land data on a rotating globe 
#  Philip Kraaijenbrink <p.d.a.kraaijenbrink@uu.nl>
#  18 Dec 2020

library(tidyverse)
library(pbapply)
library(pkrf)         # see kraaijenbrink/pkrf on github
library(raster)
library(ggspatial)
library(RStoolbox)
library(snowfall)
rm(list=ls())


# setup -------------------------------------------------------------------

# settings
outroot   = './exports/frames'
overwrite = T
flist     = list.files('./data/varstacks', full=T, patt='.tif$')
basemap   = brick('./data/basemap_4326.tiff')[[-4]]
plotext   = extent(-6430000, 6430000, -6430000, 6430000)
anisetup  = tibble(
  fn = flist,
  crmp = c('inferno','parula','tempo','magma','precip2',
           'temperature','davos','jet','wtspec','plasma',
           'turku','magma','elevnat1','spectral','globwarm',
           'viridis','bluesice','magma'),
  crmp_r = c(F,F,F,F,F,
             F,T,F,F,F,
             T,F,T,T,F,
             F,T,F)
)


# finalize setup ----------------------------------------------------------

if (dir.exists(outroot) & !overwrite){
  stop('output dir already exists')
}else if (!dir.exists(outroot)){
  dir.create(outroot)
}

# estimate variable ranges
sfInit(T, parallel::detectCores()-2, 'SOCK')
sfLibrary(raster)
sfExport('anisetup')
varranges  = sfLapply(anisetup$fn, function(x){
  cellStats(raster(x), function(y, ...){
    quantile(y, probs=c(0.02,0.98), na.rm=T)
  })
})
sfStop()
varranges     = do.call(rbind,varranges)
anisetup$low  = varranges[,1]
anisetup$high = varranges[,2]

# replicate to fit hourly data
n_var         = nrow(anisetup)
anisetup      = anisetup[rep(1:n_var, each=24),]
anisetup$band = rep(1:24, n_var)
anisetup$time = seq(as.POSIXct('2020-01-01', tz='UTC'),
                    by='hours', length.out=nrow(anisetup))
anisetup$ts   = format(anisetup$time, '%Y-%m-%d %H:%M')

# replicate for smoother rotation framerate
anisetup      = anisetup[rep(1:nrow(anisetup),each=3),]
anisetup$lon  = approx(c(1,nrow(anisetup)),c(100,-85), xout=1:nrow(anisetup))$y

# plotting ----------------------------------------------------------------

# function to generate earth outline
circleFun <- function(center=c(0,0), r=1, npoints=360){
  tt <- seq(0, 2*pi, length.out=npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(tibble(x = xx, y = yy))
}

# function to generate plot for a single frame
makePlot = function(
    fn,
    b           = 1,
    lat         = 0,
    lon         = 0,
    pext        = plotext,
    resolution  = 12000,
    colramp     = 'inferno',
    colramp_rev = F,
    l           = 240,
    h           = 310,
    bordercolor = 'white',
    timestamp   = '2020-01-01 00:00'){
  
  # reproject era data and basemap
  projdef      = sprintf('+proj=ortho +lon_0=%f +lat_0=%f', lon, lat)
  pdat         = projectRaster(raster(fn,band=b), crs=crs(projdef), res=resolution)
  pdat_pts     = rasterToPoints(pdat) %>% as_tibble %>% rename(z=3)
  
  # reproject base layer
  basemap_proj = projectRaster(basemap, crs=crs(projdef), res=resolution)
  baselayer    = ggRGB(basemap_proj,1,2,3,ggLayer=T)
  
  # clean up variable name for title
  varname = fn %>%
    str_extract('[0-9]{8}_.*.tif$') %>%
    str_remove('_hourly') %>% 
    str_remove('[0-9]{8}_') %>%
    str_remove('.tif$') %>%
    str_replace_all('_',' ') %>%
    str_to_sentence()
  
  # actual plotting
  ggplot(pdat_pts) + 
    baselayer +
    geom_raster(aes(x,y,fill=z), hjust=0.5, vjust=0.5, alpha=0.75) + 
    scale_fill_gradientn(colours=pkrf::ramp(colramp,rev=colramp_rev),
                         limits=c(l,h), oob=scales::squish, name=NULL) +
    geom_path(aes(x,y), data=circleFun(r=pext[2]), colour=bordercolor, size=2.5) +
    coord_sf(expand=F, datum=NULL, xlim=pext[1:2], ylim=pext[1:2]) + 
    theme_pk(base_size=14) + 
    labs(x=NULL, y=NULL, title=varname, subtitle=timestamp) + 
    labs(caption='ERA5-Land (ECMWF)                                                   @philipkraai') + 
    theme(plot.background      = element_blank(),
          panel.background     = element_blank(),
          panel.border         = element_blank(),
          axis.line            = element_blank(),
          plot.title           = element_text(margin=margin(b=9), hjust=0.5, face='plain', size=rel(1.3)),
          plot.subtitle        = element_text(margin=margin(b=12), hjust=0.5, face='bold', size=rel(0.9), family='mono'),
          plot.caption         = element_text(margin=margin(0,0,0,0), face='plain', size=rel(0.8)),
          legend.position      = 'none'
    )
}


# export ------------------------------------------------------------------

# export frames in parallel (RAM disk swapping may corrupt images)
invisible(gc())
sfInit(T, min(c(8,parallel::detectCores()-2)), 'SOCK')
sapply(c('raster','ggplot2','stringr','dplyr','raster','ggspatial','RStoolbox','pkrf'),
       sfLibrary, char=T)
sfExportAll()
out = sfLapply(1:nrow(anisetup), function(i){
  p = makePlot(fn=anisetup$fn[i], lat=0, lon=anisetup$lon[i], b=anisetup$band[i],
               resolution=10000, colramp=anisetup$crmp[i], colramp_rev=anisetup$crmp_r[i],
               l=anisetup$low[i], h=anisetup$high[i], timestamp=anisetup$ts[i])
  ggsave(file.path(outroot,sprintf('%04d.png',i)), p, dpi=200, width=15, height=15, units='cm')
  return(NULL)
})
sfStop()
