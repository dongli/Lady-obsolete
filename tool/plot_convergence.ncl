load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	x = (/0.375,0.75,1.5,3.0/)
	errors = (/0.00765,0.01546,0.03159,0.06551/)

	; calculate the numerical convergence rate
	log_x = log(x)
	log_e = log(errors)
	rate = regline(log_x, log_e)

	; calculate the minimal resolution
	min_e = 0.033
	log_min_e = log(min_e)
	do i = 0, dimsizes(errors)-1
		if (errors(i) .ge. min_e) then
			x_span = log_x(i)-log_x(i-1)
			e_span = log_e(i)-log_e(i-1)
			x_m = exp((log_min_e*x_span-log_x(i)*log_e(i-1)+log_x(i-1)*log_e(i))/e_span)
			break
		end if
	end do
	if (.not. isdefined("x_m")) then
		system("[Error]: Can not find out minimal resolution!")
		exit
	end if

	; reference Lines
	ref_line_1st = 0.1*x^2
	ref_line_2nd = 0.001*x^3

	wks = gsn_open_wks("pdf", "convergence")

	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@vpWidthF = 0.6
	res@vpHeightF = 0.4
	res@xyXStyle = "Log"
	res@xyYStyle = "Log"
	res@xyLineThicknessF = 3.0
	res@trXReverse = True
	res@trXMinF = min(x)-0.1
	res@trXMaxF = max(x)+0.5
	res@trYMinF = 1.0e-5
	res@trYMaxF = 1
	res@tmXBTickStartF = res@trXMinF
	res@tmXBTickEndF = res@trXMaxF
	res@tmYLTickStartF = 1.0e-8
	res@tmYLTickEndF = 1
	res@tmXMajorGrid = True
	res@tmYMajorGrid = True
	res@tmXMajorGridThicknessF = 1.0
	res@tmYMajorGridThicknessF = 1.0
	res@tmXMajorGridLineDashPattern = 2
	res@tmYMajorGridLineDashPattern = 2
	res@tmXBMode = "Explicit"
	res@tmXBValues = x
	res@tmXBLabels = (/"0.375~S~o~N~","0.75~S~o~N~","1.5~S~o~N~","3~S~o~N~"/)
	res@tiXAxisString = "resolution"
	res@tiYAxisString = "L~B~inf~N~ error"
	res@gsnLeftString = "convergence rate: "+sprintf("%4.2f", rate)

	res@xyMarkLineMode = "MarkLines"
	res@xyMarkerSizeF = 10
	plot0 = gsn_csm_xy(wks, x, errors, res)

	res@xyMarkLineMode = "Lines"
	plot1 = gsn_csm_xy(wks, x, ref_line_1st, res)
	plot2 = gsn_csm_xy(wks, x, ref_line_2nd, res)

	overlay(plot0, plot1)
	overlay(plot0, plot2)

	draw(plot0)
	frame(wks)

end