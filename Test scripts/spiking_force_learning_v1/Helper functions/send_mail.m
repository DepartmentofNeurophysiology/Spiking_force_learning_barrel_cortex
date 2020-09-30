function  send_mail( text )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

setpref('Internet','SMTP_Server','smtp-mail.outlook.com');
setpref('Internet','E_mail','hscheeres@hotmail.com');
setpref('Internet','SMTP_Username','hscheeres@hotmail.com');
setpref('Internet','SMTP_Password','');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable', 'true');
props.setProperty('mail.smtp.socketFactory.port','587');
sendmail('noahvandebunt@outlook.com ',text) ;

end

