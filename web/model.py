import uuid

from peewee import SqliteDatabase, Model, UUIDField, BooleanField, TextField

db = SqliteDatabase('data.db')


class Job(Model):
    class Meta:
        database = db

    id = UUIDField(primary_key=True, default=uuid.uuid4())
    fasta = TextField()
    algorithm = TextField()
    done = BooleanField(default=False)


def prepare_db():
    db.connect()
    db.create_table(Job)


if __name__ == '__main__':
    prepare_db()
