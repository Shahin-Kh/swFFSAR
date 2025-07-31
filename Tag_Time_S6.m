function time1=Tag_Time_S6(time)
% Using this function we calculate time
% in UTC Georgian format from dataset time tags

epoch = datetime(2000,01,01);
time1 = epoch + seconds(time);
time1 = datestr(time1, 'DD-mmm-YYYY');

end