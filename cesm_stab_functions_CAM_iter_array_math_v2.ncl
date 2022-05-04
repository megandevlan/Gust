load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;
; RBN: 13Dec2017: Reveal structure of surface flux stability relationships.
;      20Dec2017: Read CAM 3hourly instantaneous data and plot info of 
;                 the # of iteractions required for different ocean regions.  

;;; Ideas ;;;
;; 
;   1. hol (H/L) is limited to 10, could increase?
;   2. Could make Cd a stronger function u10m?
;   3. When delt is very large 



begin


; Data to read in from run.
; Data available at 3hrly instantaneous(U,V,TS,T)

  ;run_name = "F2010climo.f09_f09_mg17.TestGustGlobal_Control_32Lcam6_branch.003"
  run_name = "F2010climo.f09_f09_mg17.TestGustGlobal_Control_58Lcam6_branch.003"

  dir_in = "/glade/scratch/mdfowler/archive/"+run_name+"/atm/hist/"

  ; it0 = 1   ; Time samples.
  ; it1 = 55

  conv_var = "tau"  ;lh/sh/tau/hol

  fileOut = "/glade/work/mdfowler/gust_cam7/Analysis/"+run_name+".stabilityInfo.nc"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
   fileName = dir_in+run_name+".cam.h1.0004-01-01-10800.nc"
   fileIn   = addfile(fileName, "r") 

  tbot_loop = fileIn->TBOT(:,:,:)
  ts_loop   = fileIn->TS(:,:,:)
  ubot_loop = fileIn->UBOT(:,:,:)
  vbot_loop = fileIn->VBOT(:,:,:)
  qbot      = fileIn->QBOT(:,:,:) 
  ; zbot      = fileIn->ZBOT(:,:,:)
  thbot     = fileIn->THBOT(:,:,:)
  rbot      = fileIn->RBOT(:,:,:)  

  lat = fileIn->lat
  lon = fileIn->lon
  time = fileIn->time
  date = fileIn->date
  ; lev  = fileIn->lev(it0:it1)

  nlats = dimsizes(lat)
  nlons = dimsizes(lon)
  ntime = dimsizes(time)
  ; nlev  = dimsizes(lev)
 
  fmiss = default_fillvalue(typeof(tbot_loop))



; Reading in OCEFRAC
  file_ofrac = "/glade/work/mdfowler/gust_cam7/output/F2010climo.f09_f09_mg17.TestGustGlobal_Control.001.cam.h0.0002-01.nc"
  fin_ofrac = addfile(file_ofrac,"r")
  ofrac = fin_ofrac->OCNFRAC(0,:,:)

;; Vars copied from CAM just to be sure
;
;  SHR_CONST_BOLTZ   = 1.38065e-23
;  SHR_CONST_AVOGAD  = 6.02214e26 
;  SHR_CONST_MWWV    = 18.016  
;  SHR_CONST_MWDAIR  = 28.966
;  SHR_CONST_RGAS    = 6.02214e26*1.38065e-23;
  
  SHR_CONST_RDAIR   = (6.02214e26*1.38065e-23)/28.966 
  SHR_CONST_RWV     = (6.02214e26*1.38065e-23)/18.016

  loc_zvir = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1

  loc_karman = 0.4
  loc_cpvir = 1.81e3/1.00464e3
  loc_cpdair = 1.00464e3

  loc_g = 9.80616
  loc_latvap = 2.501e6 ; LH f evap.
  loc_stebol = 5.67e-8 ; SB constant.

  zref = 10. ; Ref height for U
  ztref = 2. ; Ref height for T
 
  ; !!!!!!! 
  ; These two values should change with resolution: 
  ; zbot = 66.
  ; pbot = 992.556 ; PBOT (mb)
  zbot = 22. 
  pbot = 997.528

  ; zbot = 50.  ; Lowest model level.
  umin = 0.5 ; Floor for vmag.
  ocn_frac_lim = 0.01 ; Ofrac has to be greater than this for calcs.
  
  al2 = log(zref/ztref)

  vmagx_loop = 0 ; Standard vmag used in calculation.
  ; vmagx_loop = 1 ; Standard vmag used in calculation.
  hol_lim0 = 10. ; Stab and instab limiter on hol
 
  niter = 5   ; # of iterations in default shr_flux_mod, "flux_con_max_iter"
  ; niter = 20 ;# of implicit iterations.
  thresh_pcent = 0.1  ; Convergence criteria %.

  pnon_conv = 0 ; Non-converging point counter


;;;; Set up arrays ;;;;
  
  nconv_iter = 0*tointeger(tbot_loop)
  imiss =  default_fillvalue(typeof(nconv_iter))

  nconv_iter = default_fillvalue(typeof(nconv_iter))
  
  flux_conv = new((/niter+1,ntime,nlats,nlons/),"float")

;;;;;;;;;;;;;;; ARRAY LOOPS ;;;;;;;;;;;;;;;;;

  nloops = dimsizes(ts_loop)

;  do itime=0,ntime-1
;    print("itime = "+(itime+1)+" of "+ntime)
;    do ilat=0,nlats-1
;      do ilon=0,nlons-1
        
;        if (ofrac(ilat,ilon).gt.ocn_frac_lim) then
        
; Assign array values to this loop
        
          
;        qbot = qbot_loop
          ts = ts_loop  ; TS (SST)          
          tbot = tbot_loop ; TBOT
          ubot = ubot_loop ; UBOT
          vbot = vbot_loop ; VBOT
          
          print("-Assign up front info")

; Derived.
          ; thbot = tbot*(pbot/1000.)^0.286 ;THBOT (pot. temp) form from ncl.
          ; rbot = pbot*100./(8.314e3*tbot)      ; Air density - pbot/(Rsp*T_bot)
          
          
; Derived state settings
         
          qsat = 640380.0 / exp(5107.4/ts)
          ;   MDF: Read in qbot above, shouldn't need this approximation now 
          ; qbot = 0.90*qsat  ; RBN: Approximation given that we don't have qbot data.
          thvbot = thbot * (1.0 + loc_zvir * qbot) ; virtual potential temp (K)

          ssq = 0.98 * qsat / rbot
          delt = thbot - ts 
          delq = qbot - ssq 
          alz    = log(zbot/zref) 
          cp     = loc_cpdair*(1.0 + loc_cpvir*ssq) 
          hol_limit = hol_lim0+ts*0.
          
;;;  Mahrt and Sun (1995) ;;;
; Cold Air Outbreak Modification: 
; Increase windspeed for negative tbot-ts 

          alpha = 1.4 ; 
          maxscl = 2.0 ; maximum wind scaling for flux
          td0 = -10.  ; start t-ts for scaling
          
;          vmag0 = max((/umin, sqrt(ubot^2+vbot^2)/))
          vmag0 = where(sqrt(ubot^2+vbot^2).gt.umin,sqrt(ubot^2+vbot^2),umin)

         ; if (delt.lt.td0) then 
            ;vscl  = min((/(1.+alpha*(abs(delt-td0)^0.5/abs(vmag0))),maxscl/))
         ;    vscl  = where((1.+alpha*(abs(delt-td0)^0.5/abs(vmag0))).lt.maxscl,(1.+alpha*(abs(delt-td0)^0.5/abs(vmag0))),maxscl)
         ; else
         ;   vscl = 1.
         ; end if
          
;; Options to apply to wind.
          
          if (vmagx_loop.eq.0) then ; Standard
            vx = 1.
            vy = 1.
          end if
          
          if (vmagx_loop.eq.1) then ; Apply extra wind speed (vmag) for alll calc.
            vx = vscl
            vy = vscl
          end if
          
          if (vmagx_loop.eq.2) then ; Apply extra wind speed to endd flux only
            vx = 1. 
            vy = vscl
          end if
          
          vmag = vmag0*vx
          
          
; Z/L ustar, tstar first estimates
          
          stable = 0.5 + sign_f90(delt*0.+0.5 , delt) ; Basically routines 0.5 with the sign of delt.
          
          cdn    = 0.0027 / vmag + 0.000142 + 0.0000764 * vmag ; Neutral drag coeff at lowest model layer wind value.
          rdn    = sqrt(cdn)  ; sqrtt of neutral exchange coeff (momomentum)
          
; Neutral exchange coeff (iff stable=0.0327, unstable=0.018)
          rhn    = (1.0-stable) * 0.0327 + stable * 0.018 ; sqrtt of neutral exchange coeff (heat)
          ren    = 0.0346 ; Water.
          
; Friction quenstities (shear velocity etc.)
          ustar = rdn * vmag
          tstar = rhn * delt  
          qstar = ren * delq  
          
; State output.
          
;          print("")
;          print("####    Location (lat="+sprintf("%5.2f",lat(ilat))+", lon="+sprintf("%5.2f",lon(ilon))+")"+": TS="+ts+", TBOT="+thbot+", DT="+(ts-thbot)+", VMAG="+vmag+", QBOT="+(1000.*qbot)+", QSAT="+(1000.*qsat)+", HOLlim="+hol_limit+", VSCL="+vscl+" (vx="+vx+",vy="+vy+")   ####")

   




;;;;;;
;;;; Compute stability and evaluate alll stability functions iterative solutions.    
;;;;;;   

          do iter=0,niter
            print("iter = "+(iter+1)+" of "+niter)


; H/L at Zbot "stability parameter"
            hol  = loc_karman*loc_g*zbot*(tstar/thbot+qstar/(1.0/loc_zvir+qbot))/ustar^2 ; (H/L at Zbot)
;      print("hol="+hol)


            hol  = sign_f90(dim_min_n((/abs(hol),hol_limit/),0), hol ) ; Limits hol (H/L) to 10x

; Determine whether this is stable or unstable.
            stable = 0.5 + sign_f90(0.5+0.*hol , hol) ; Iff stable = 1, (+hov), stable: iff stable = 0 (-hov) unstable. 

            xsq    = dim_max_n((/sqrt(abs(1.0+0.*hol - 16.0*hol)) , 1.0+0.*hol/),0) ; Factor for the flux profile function
            xqq    = sqrt(xsq) ; ""

; Flux profiles (iff unstable)
            psimhu = log((1.0+xqq*(2.0+xqq))*(1.0+xqq*xqq)/8.0) - 2.0*atan(xqq) + 1.571
            psixhu = 2.0 * log((1.0 + xqq*xqq)/2.0)
            
; Stability function ats at zbot 
; Selection of stability parameter based on hol (iff stable = 1, stable (first term): iff stable = 0, unstable (second term).
            psimh  = -5.0*hol*stable + (1.0-stable)*psimhu
            psixh  = -5.0*hol*stable + (1.0-stable)*psixhu

; Shift wind speed using old coefficient
            rd   = rdn / (1.0 + rdn/loc_karman*(alz-psimh))
            u10n = vmag * rd / rdn 

; Update transfer coeffs at 10m and neutral stability 
            cdn = 0.0027 / u10n + 0.000142 + 0.0000764 * u10n ; Neutral drag coeff at 10m.
            rdn = sqrt(cdn)
            ren = 0.0346
            rhn = (1.0-stable)*0.0327 + stable * 0.018 
            
; Shift alll coeffs to measurement height and stability 
            rd = rdn / (1.0 + rdn/loc_karman*(alz-psimh)) 
            rh = rhn / (1.0 + rhn/loc_karman*(alz-psixh)) 
            re = ren / (1.0 + ren/loc_karman*(alz-psixh)) 
;      print("++ PSIMH="+psimh+", PSIXH="+psixh+", ln(Za/Zo)="+alz)
          
; Update ustar, tstar, qstar using updated, shifted coeffs

            vmag = vmag0*vy ; Iff add component to flux calculation only.
            
            ustar = rd * vmag 
            tstar = rh * delt 
            qstar = re * delq 
            
; Compute the final fluxes
            
            tau = rbot * ustar * ustar 
            
; Momentum flux ---
            taux = tau * ubot / vmag 
            tauy = tau * vbot / vmag  
            
; Heat flux ---
            sen_hf  =          cp * tau * tstar / ustar 
            lat_hf  =  loc_latvap * tau * qstar / ustar
            lwup = -loc_stebol * ts^4 
            
; Water flux ---
            evap = lat/loc_latvap 
            
            vmag = vmag0*vx
          
            ustar = rd * vmag 
            tstar = rh * delt 
            qstar = re * delq 
          
          
;            print("iter = "+(iter+1)+" ** Stable="+stable+", hol="+hol+", SH="+sprintf("%5.2f",sen_hf)+", LH="+sprintf("%5.2f",lat_hf)+", TAUX="+sprintf("%5.3E",taux)+", TAUY="+sprintf("%5.3E",tauy)+", U10="+sprintf("%5.2f",u10n))     
          
; Save variables to be (maybe) plotted

            if (conv_var.eq."lh") then
              flux_conv(iter,:,:,:) = lat_hf
            else if (conv_var.eq."sh") then
              flux_conv(iter,:,:,:) = sen_hf
            else if (conv_var.eq."tau") then
              flux_conv(iter,:,:,:)    = tau
            end if
            end if
            end if
          

; Update convergence iteration #
 ;           if (iter.gt.0) then
 ;             if (100.*(flux_lat(iter)-flux_lat(iter-1))/flux_lat(iter-1).lt.thresh_iter) then
 ;               if (ismissing( nconv_iter(itime,ilat,ilon))) then ; Apply convergence iteration # (iff not already reached).
 ;                 nconv_iter(itime,ilat,ilon) = iter
 ;               end if
 ;             end if
 ;           end if
            

          end do ; ITERATION

; Iff does not converge.
;          if (ismissing(nconv_iter(itime,ilat,ilon))) then
;            print("Non converging point at : "+date(itime)+", lat="+lat(ilat)+", lon="+lon(ilon))
;            pnon_conv = pnon_conv+1
;          end if

; Diagnotic calculations for 2m
          hol = hol*ztref/zbot
          xsq = dim_max_n( (/1.0+0.*hol, sqrt(abs(1.0+0.*hol-16.0*hol))/),0)
          xqq = sqrt(xsq)
          psixhu = 2.0 * log((1.0 + xqq*xqq)/2.0)
          psix2   = -5.0*hol*stable + (1.0-stable)*psixhu
          fac     = (rh/loc_karman) * (alz + al2 - psixh + psix2 )
          tref = thbot - delt*fac 
          qref = qbot - delq*fac
          
;          print("T2m="+tref+", Q2m="+(1000.*qref))
          
;        else
; Not an ocean point
;         print("####    Location (lat="+sprintf("%5.2f",lat(ilat))+", lon="+lon(ilon)+"): Non-ocean gird point") 
;        end if
 
;      end do ; TIME
;    end do ; LAT
;  end do ; LON
  

;;;; Post calculations
;;
;; Write out some info to .nc file 

; ----- Create netCDF file ------
print("Creating netCDF file")

system("/bin/rm -f " + fileOut)   ; remove any pre-existing file
ncdf = addfile(fileOut,"c")  ; open output netCDF file

; Create global attributes of the file (optional)
; fAtt   = False
fAtt               = True            ; assign file attributes
fAtt@title         = "Offline calculation of variables in shr_flux_mod"
fAtt@source_file   =  fileName
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

; Make time an unlimited dimension
filedimdef(ncdf,"time",-1,True)

;===================================================================
; output variables directly; NCL will call appropriate functions
; to write the meta data associated with each variable
;===================================================================
u10n!0 = "time"
u10n!1 = "lat"
u10n!2 = "lon"

u10n&time  = fileIn->time
u10n&lat  = fileIn->lat
u10n&lon  = fileIn->lon

u10n@long_name = "u10n value in shr_flux_mod"
u10n@units = "m/s"
ncdf->u10n = u10n


tref!0 = "time"
tref!1 = "lat"
tref!2 = "lon"

tref&time  = fileIn->time
tref&lat  = fileIn->lat
tref&lon  = fileIn->lon

tref@long_name = "tref value in shr_flux_mod"
tref@units = "K"
ncdf->tref = tref



stable!0 = "time"
stable!1 = "lat"
stable!2 = "lon"

stable&time  = fileIn->time
stable&lat  = fileIn->lat
stable&lon  = fileIn->lon

stable@long_name = "Stability value in shr_flux_mod"
stable@units = "-"
ncdf->stable = stable


hol!0 = "time"
hol!1 = "lat"
hol!2 = "lon"

hol&time  = fileIn->time
hol&lat  = fileIn->lat
hol&lon  = fileIn->lon

hol@long_name = "Value of H/L in shr_flux_mod"
hol@units = "-"
ncdf->hol = hol



psimh!0 = "time"
psimh!1 = "lat"
psimh!2 = "lon"

psimh&time  = fileIn->time
psimh&lat  = fileIn->lat
psimh&lon  = fileIn->lon

psimh@long_name = "PSIMH value in shr_flux_mod"
psimh@units = "-"
ncdf->psimh = psimh



rd!0 = "time"
rd!1 = "lat"
rd!2 = "lon"

rdtime  = fileIn->time
rd&lat  = fileIn->lat
rd&lon  = fileIn->lon

rd@long_name = "rd value in shr_flux_mod"
rd@units = "-"
ncdf->rd = rd


rdn!0 = "time"
rdn!1 = "lat"
rdn!2 = "lon"

rdntime  = fileIn->time
rdn&lat  = fileIn->lat
rdn&lon  = fileIn->lon

rdn@long_name = "rdn value in shr_flux_mod"
rdn@units = "-"
ncdf->rdn = rdn

;; Iteration to convergence

  print("Number of non-converging points = "+pnon_conv+" out of "+(ntime*nlats*nlons)+" total points")
  
;  nconv_iter = tointeger(tbot_loop)
;  printVarSummary(nconv_iter)

  print(ind(ndtooned(flux_conv.eq.0.)))

  do iter=1,niter
    flux_diter = 100.*abs((flux_conv(iter,:,:,:)-flux_conv(iter-1,:,:,:))/flux_conv(iter,:,:,:)) ; Convergence criteria met?
    nconv_iter = where(flux_diter.lt.thresh_pcent.and.nconv_iter.eq.imiss,iter,nconv_iter) ; Set iteration number for this point.
  end do
  

;  nconv_iter = where(nconv_iter.eq.imiss,-100,nconv_iter) ; Set to missing value for no convergence

  nconv_iter@_FillValue = imiss


;; Average/maxx iteration #

  nconv_iter_ave = dim_avg_n(nconv_iter,0) 
  nconv_iter_max = dim_max_n(nconv_iter,0) 

  
  copy_VarMeta(ofrac,nconv_iter_ave) 
  copy_VarMeta(ofrac,nconv_iter_max) 

; %age of non converging points.

  nconv_miss = tofloat(nconv_iter_ave)

  do ilat=0,nlats-1
    do ilon=0,nlons-1
      ii := ind(ismissing(nconv_iter(:,ilat,ilon)))
      nconv_miss(ilat,ilon) = tofloat(dimsizes(ii))
      if (ismissing(nconv_miss(ilat,ilon))) then
        nconv_miss(ilat,ilon) = 0. ; Set to 0 if ind=imiss value
      end if
    end do
  end do



 
  nconv_miss = 100.*nconv_miss/tofloat(ntime) ; %age of total.
;  print(nconv_miss)
  printMinMax(nconv_miss,True)

  copy_VarMeta(ofrac,nconv_miss) 


end 
