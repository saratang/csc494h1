classdef fifo < handle % why do I have to implement this
    properties
        head
        tail % a pointer (?) to the tail for easy appending
    end
    methods
        function obj = fifo()
            obj.head = node();
            obj.tail = obj.head;
        end
        
        function add(obj, val)
            newNode = node(val);
            obj.tail.next = newNode;
            obj.tail = obj.tail.next;
        end
        
        function n = peek(obj)
            n = obj.head.val;
        end
        
        function n = pop(obj)
            n = obj.head.val;
            obj.head = obj.head.next;
        end
        
        function e = empty(obj)
            e = ISNAN(obj.head.val);
        end
    end
end