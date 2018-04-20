function dbConn = connect2TotalsDB
    %Connect to MySQL coolops database
    host = 'localhost';
    user = 'root';
    password = 'password';
    dbName = 'dBname';
    jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
    jdbcDriver = 'com.mysql.jdbc.Driver';
    % Make the connection
    dbConn = database(dbName, user , password, jdbcDriver, jdbcString);
end