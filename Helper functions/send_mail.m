function  send_mail( text )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','timwinter06@gmail.com');
setpref('Internet','SMTP_Username','timwinter06');
setpref('Internet','SMTP_Password','timber1995');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('timwinter06@gmail.com',text) ;

end

