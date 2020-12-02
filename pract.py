class MyClass:


    def __init__(self, name, number):
        self.name = name
        self.number = number

    def __str__(self):
        return f"the object= {self.name}, {self.number}"

    def __repr__(self):
        return f"<the object= {self.name}, {self.number}>"

MyClass("name1", 1)
print(MyClass)
a = MyClass("name1", 1)
print("a=", a)
b = MyClass("name1", 2)
print("a=", a)
print("b=", b)
