load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	x = fspan(1, 4, 4)
	errors = (/0.00262,0.00524,0.00717,0.01071/)

	; reference Lines
	ref_line_1st = 0.01*x
	ref_line_2nd = 0.0001*x^2

	wks = gsn_open_wks("pdf", "convergence")

	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@vpWidthF = 0.4
	res@vpHeightF = 0.3
	res@xyXStyle = "Log"
	res@xyYStyle = "Log"
	res@xyLineThicknessF = 5.0
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
	res@tiXAxisString = "temporal resolution"
	res@tiYAxisString = "Linf error"
	res@gsnLeftString = "convergence rate: "+sprintf("%6.4f", rate)

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