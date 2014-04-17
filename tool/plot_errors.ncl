load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    start_time = 0
    end_time = 576
    time_step = 1
    file_prefix1 = "tracers.rotation.eq.true.240x120"
    file_prefix2 = "tracers.rotation.eq.240x120"

    fs1 = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/Release/"+file_prefix1+".%3.3d.nc\n' $i; done")

    fs2 = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/Release/"+file_prefix2+".%3.3d.nc\n' $i; done")

    if (dimsizes(fs1) .ne. dimsizes(fs2)) then
        system("echo [Error]: The sizes of files are not matched!")
    end if
    num_time = dimsizes(fs1)

    errors = new((/3,num_time/), "double")
    do l = 0, num_time-1
        system("echo file: "+fs1(l))
        f1 = addfile(fs1(l), "r")
        f2 = addfile(fs2(l), "r")

        dq = abs(f1->q0-f2->q0)
        errors(0,l) = sum(dq*f1->volume)/sum(abs(f1->q0)*f1->volume)
        errors(1,l) = sqrt(sum(dq^2*f1->volume)/sum(f1->q0^2*f1->volume))
        errors(2,l) = max(dq*f1->volume)/max(abs(f1->q0)*f1->volume)

        system("echo "+sprintf("%20.5f", errors(0,l))+sprintf("%20.5f", errors(1,l))+sprintf("%20.5f", errors(2,l)))

    end do

    wks = gsn_open_wks("pdf", file_prefix2+".error")

    res = True
    res@xyLineThicknessF = 3
    res@xyDashPatterns = (/0,1,2/)
    res@trXMinF = 0
    res@trXMaxF = 576
    res@trYMinF = 0
    res@trYMaxF = 0.1
    res@tmXBTickStartF = 0
    res@tmXBTickEndF = 576
    res@tmYLTickStartF = 0
    res@tmYLTickEndF = 0.1

    plot = gsn_csm_xy(wks, ispan(0, num_time-1, 1), errors, res)

end