def twoSum(number_list,target):
    '''
    :param number_list: list of numbers
    :param target: target sum
    :return: index of two numbers that sum up to target
    '''
    if target % 2 == 0:
        half_target = target/2
        indexes = []
        for k in range(len(number_list)):
            if number_list[k] == half_target:
                indexes.append(k)
        if len(indexes) == 2:
            return indexes
    other_val = [target-i for i in number_list]
    numbers = set(other_val).intersection(set(number_list))
    indexes = []
    for each_val in numbers:
        if each_val != (target/2):
            indexes.append(number_list.index(each_val))
    return(indexes)

def display(node):
    while node is not None:
        print(node.val)
        node = node.next

class ListNode:
     def __init__(self, val=0, next=None):
         self.val = val
         self.next = next

class Solution:
    def addTwoNumbers(self, list1, list2):
        list_1 = []
        while list1 is not None:
            list_1.append(list1.val)
            list1 = list1.next
        number_1 = "".join(map(str,list_1))[::-1]
        list_2 = []
        while list2 is not None:
            list_2.append(list2.val)
            list2 = list2.next
        number_2 = "".join(map(str,list_2))[::-1]
        sum = str(int(number_1) + int(number_2))
        result = list(map(int,sum))[::-1]
        node = new_node = ListNode(result[0])
        for i in result[1:]:
            temp = ListNode(i)
            new_node.next = temp
            new_node = new_node.next
        return(node)

    def lengthOfLongestSubstring(self,str):
        if len(str) == 0:
            return 0
        def countlength(str):
            count = 0
            s = set()
            for i in str:
                if i not in s:
                    count += 1
                    s.add(i)
                else:
                    break
            return count
        length_count = []
        for i in range(len(str)):
            length_count.append(countlength(str[i:]))
        return(max(length_count))

    def findMedianSortedArrays(self, nums1, nums2):
        '''
        run time should be log(m+n)
        :param nums1:
        :param nums2:
        :return: float
        '''
        total = len(nums1) + len(nums2)
        def getSmaller(list1,list2):
            if len(list1) == 0:
                list2 = list2[1:]
            elif len(list2) == 0:
                list1 = list1[1:]
            elif list1[0]<=list2[0]:
                list1 = list1[1:]
            else:
                list2 = list2[1:]
            return(list1,list2)
        count = 0
        while count < total/2-1:
            nums1, nums2 = getSmaller(nums1,nums2)
            count +=1
        if (total) % 2 ==1:
            try:
                return(float(min(nums1[0],nums2[0])))
            except:
                if len(nums1)==0:
                    return nums2[0]
                else:
                    return nums1[0]
        else:
            try:
                first = min(nums1[0],nums2[0])
            except:
                if len(nums1)==0:
                    first = nums2[0]
                else:
                    first = nums1[0]
            nums1, nums2 = getSmaller(nums1,nums2)
            try:
                second = min(nums1[0], nums2[0])
            except:
                if len(nums1) == 0:
                    second = nums2[0]
                else:
                    second = nums1[0]
            return(float((first+second)/2))

    def longestPalindromicSubstring(self,str):
        def longestPalindrome(self, str):
            for i in range(len(str), 0, -1):
                for k in range(0, int(len(str) - i + 1)):
                    candidate = str[k:k + i]
                    if candidate == candidate[::-1]:
                        return candidate

    def convertZigZag(self,str,nrow):
        if nrow == 1:
            return(str)
        content = {}
        #initate row content
        for i in range(nrow):
            content[i]=""
        print(content)
        count = 0
        i = -1
        for index, char in enumerate(str):
            content[count] += char
            if count == (nrow-1) or count == 0:
                i = -i
            count += i
        return("".join(content.values()))

    def


sol = Solution()
print(sol.convertZigZag("AB",1))








