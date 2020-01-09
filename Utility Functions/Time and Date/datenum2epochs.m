function epochs = datenum2epochs(datenums)



[Y,M,D,H,MN,S] = datevec(datenums);


epochs = jd2epochs(cal2jd(Y,M,D+H/24+MN/(60*24)+S/86400));

end
