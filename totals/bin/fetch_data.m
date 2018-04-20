function result = fetch_data(session, statement)
% Query database for data on system type
if isopen(session)
    result = get(fetch(exec(session, statement)), 'Data');
else
    fprintf(1, 'Connection failed: %s', session.Message);
    return
end