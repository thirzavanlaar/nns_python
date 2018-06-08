;---------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;---------------------------------------------------------------
begin

File = addfile("subsubgrid.nc","r")
File2 = addfile("subsubdomain_lasttimestep.nc","r")
;File = addfile("/data/inscape/vanlaar/HDCP2_data/input/TA_dom4_qc_time41.nc","r")
ncells = getfilevardimsizes(File,"clon");

; read grid information (i.e. coordinates of cell centers and vertices)
rad2deg = 45./atan(1.)      ; radians to degrees
clon = File->clon * rad2deg ; cell center, lon (ncells)
clat = File->clat * rad2deg ; cell center, lat (ncells)

vlon = File2->clon_bnds * rad2deg   ; cell vertices, lon (ncells,3)
vlat = File2->clat_bnds * rad2deg   ; cell vertices, lat (ncells,3)

x = clon
y = clat

; read qc data

inputfile = addfile("subsub_regular.nc","r")

qc = File2->qc(:,:,:)   ; dims: (time,lev,cell)

qc_regular1 = inputfile->qc(0,0,:)

qc_regular = floattointeger(qc_regular1)

print(sum(qc_regular))


qc_mask = new(ncells,"integer")
qc_mask(:) = 0

do i=0,ncells-1
  column = sum(qc(0,:,i))
  if (column.gt.0.) then
    qc_mask(i) = 1
  end if
end do





;-------------------------------------------------------------------------
; create plot

wks = gsn_open_wks("pdf","../Figures/cloudmask_regular_ncl")
gsn_define_colormap(wks,"testcmap")      ; choose colormap

ResC                     = True
ResC@gsnFrame	         = False
ResC@gsnDraw	         = False
ResC@gsnMaximize         = True
ResC@cnInfoLabelOn       = False
ResC@lbLabelBarOn	 = False

ResC@sfXArray            = clon        ; cell center (lon)
ResC@sfYArray            = clat        ; cell center (lat)
ResC@sfXCellBounds       = vlon        ; define triangulation
ResC@sfYCellBounds       = vlat        ; define triangulation

ResC@cnFillOn            = True        ; do color fill
ResC@cnFillMode          = "cellfill"
ResC@cnLinesOn           = False       ; no contour lines
;ResC@tiYAxisString	  = "" ;"latitude"
;ResC@tiXAxisString	  = "" ; "longitude"

ResC@pmTickMarkDisplayMode = "Always"


ResC@cnFillColors = (/"White","lightskyblue"/)
ResC@cnLevelSelectionMode = "Explicitlevels"
ResC@cnLevels = (/0.5/)

ResC@cnLineLabelsOn = False
ResC@cnLinesOn = False

ResC@tmXTOn              = False
ResC@tmYROn              = False
ResC@tmXBMinorOn         = False
ResC@tmYLMinorOn         = False

ResC@mpMinLonF           = -57.
ResC@mpMaxLonF           = -56.5
ResC@mpMinLatF           = 13.3
ResC@mpMaxLatF           = 13.5
ResC@mpFillOn = False

restext = True
restext@txFontHeightF = 0.002



plot = gsn_csm_contour_map(wks,qc_regular,ResC)
idx = gsn_add_text(wks,plot,"0",clon(0),clat(0),restext)
idx = gsn_add_text(wks,plot,"1",clon(1),clat(1),restext)
idx = gsn_add_text(wks,plot,"2",clon(2),clat(2),restext)
idx = gsn_add_text(wks,plot,"3",clon(3),clat(3),restext)
idx = gsn_add_text(wks,plot,"4",clon(4),clat(4),restext)
idx = gsn_add_text(wks,plot,"5",clon(5),clat(5),restext)
idx = gsn_add_text(wks,plot,"6",clon(6),clat(6),restext)
idx = gsn_add_text(wks,plot,"7",clon(7),clat(7),restext)
idx = gsn_add_text(wks,plot,"8",clon(8),clat(8),restext)
idx = gsn_add_text(wks,plot,"9",clon(9),clat(9),restext)


draw(plot)
frame(wks)


end

