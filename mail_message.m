function mail_message(word)
    % your account and password.
    mail_gmail = 'naoya.inoue.6.28.1024.eagle@gmail.com';
    passward = 'clcvdtqhucxgnvsy';
    % setting gmail account
    setpref('Internet', 'E_mail', mail_gmail);
    setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
    setpref('Internet', 'SMTP_Username', mail_gmail);
    setpref('Internet', 'SMTP_Password', passward); % アプリパスワード
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');

    
    target_mail = 'naoya.inoue.6.28.1024.eagle@gmail.com';
    cl_fin = clock;
    sendmail(target_mail,...
            'Matlab Message',...
            strcat('Sand mail time at ',...
                    num2str(cl_fin(1)),'-',num2str(cl_fin(2)),'-',num2str(cl_fin(3)),'-',...
                    num2str(cl_fin(4)),'-',num2str(cl_fin(5)),'-',num2str(cl_fin(6)),...
                    newline,...
                    word...
                )...
            )
end

