from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

from .models import metadata


class DB():
    def __init__(self, db_file):
        self.url = 'sqlite:///' + db_file
        self.engine = None
        self.session = None

    def conn(self):
        try:
            self.engine = create_engine(self.url, echo=False)
            Session = sessionmaker(self.engine)
            self.session = Session()
        except Exception as e:
            print(str(e))
            return False

        try:
            metadata.create_all(self.engine)
        except Exception as e:
            print(str(e))
            return False
        else:
            return True

    def close(self):
        self.session.close()
        self.engine.dispose()
