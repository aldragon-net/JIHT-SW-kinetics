class McKennaBurner:
    def __init__(self, T_room, T_burner, T_body, height):
        self.T_room = T_room
        self.T_burner = T_burner
        self.T_body = T_body
        self.height = height


burner = McKennaBurner(
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600
)
