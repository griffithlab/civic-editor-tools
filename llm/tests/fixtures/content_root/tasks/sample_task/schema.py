from pydantic import BaseModel


class SampleOutput(BaseModel):
    answer: str


SCHEMA = SampleOutput
