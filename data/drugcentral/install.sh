sudo su - postgres -c 'createdb -T template0 drugcentral'
sudo su - postgres -c "gunzip -c $(pwd)/drugcentral.dump.08262018.sql.gz | psql drugcentral"
sudo su - postgres -c 'createuser drugcentral -w'
sudo -u postgres psql
grant all privileges on database drugcentral to drugcentral;
alter user drugcentral with encrypted password password;
\c drugcentral
GRANT USAGE ON SCHEMA public TO drugcentral;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO drugcentral;
