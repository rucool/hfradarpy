function dbConn = configs_database
    %Connect to MySQL cool_ais database
    host = 'mysql1.marine.rutgers.edu';
    user = 'coolops';
    password = 'SjR9rM';
    dbName = 'coolops';
    jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
    jdbcDriver = 'com.mysql.jdbc.Driver';
    % Make the connection
    dbConn = database(dbName, user , password, jdbcDriver, jdbcString);
end